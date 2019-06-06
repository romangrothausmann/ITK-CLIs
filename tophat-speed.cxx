////program to mimic tophat-speed.r from SITK-CLIs in order to be less mem hungry
//01: based on tophat.cxx


#include <complex>

#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkWhiteTopHatImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkImageFileWriter.h>



template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){

    typedef InputPixelType  OutputPixelType;

    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;

    int CompChunk= atoi(argv[3]);
    bool noSDI= CompChunk <= 1; // SDI only if CompChunk > 1

    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(argv[1]);
    // reader->ReleaseDataFlagOn(); // needed by filter and th
    if(noSDI){
	FilterWatcher watcherI(reader);
	watcherI.QuietOn();
	watcherI.ReportTimeOn();
	try{
	    reader->Update();
	    }
	catch(itk::ExceptionObject &ex){
	    std::cerr << ex << std::endl;
	    return EXIT_FAILURE;
	    }
	}
    else{
	reader->UpdateOutputInformation();
	}
	
    const typename InputImageType::Pointer& input= reader->GetOutput();

    double bglevel= atof(argv[4]);
    double sigma= atof(argv[5]);

    typedef itk::BinaryBallStructuringElement<InputPixelType, Dimension> StructuringElementType;
    StructuringElementType  structuringElement;

    structuringElement.SetRadius(1); // SITK procedural interface uses a default radius of 1: https://itk.org/SimpleITKDoxygen/html/namespaceitk_1_1simple.html#a106c90567c31d49b3d08fff7be93483f
    structuringElement.CreateStructuringElement();

    typedef itk::WhiteTopHatImageFilter<InputImageType, OutputImageType, StructuringElementType> FilterType;
    typename FilterType::Pointer filter= FilterType::New();
    filter->SetInput(input);
    filter->SetKernel(structuringElement);
    // filter->ReleaseDataFlagOn(); // needed by stat and add

    if(noSDI){
	FilterWatcher watcher1(filter);
	try{
	    filter->Update();
	    }
	catch(itk::ExceptionObject &ex){
	    std::cerr << ex << std::endl;
	    return EXIT_FAILURE;
	    }
	}

    const typename OutputImageType::Pointer& wth= filter->GetOutput();

    typedef itk::StatisticsImageFilter<InputImageType> StatType;
    typename StatType::Pointer stat= StatType::New();
    stat->SetInput(wth);
    if(!noSDI)
	stat->SetNumberOfStreamDivisions(CompChunk);
    FilterWatcher watcherS(stat);
    stat->Update();
    
    double MX1= stat->GetMaximum();
    double bgconst= MX1 * bglevel;
    
    typedef itk::BinaryThresholdImageFilter<InputImageType, InputImageType> ThreshType;
    typename ThreshType::Pointer th= ThreshType::New();
    th->SetInput(input); // needs re-exec reader->Update() if reader->ReleaseDataFlagOn() https://insightsoftwareconsortium.atlassian.net/browse/ITK-3351?attachmentOrder=desc
    th->SetUpperThreshold(0);
    th->SetOutsideValue(0);
    th->SetInsideValue(bgconst);
    // th->ReleaseDataFlagOn();
    FilterWatcher watcherT(th);

    typedef itk::AddImageFilter<InputImageType, InputImageType, InputImageType> AddType;
    typename AddType::Pointer add = AddType::New();
    add->SetInput1(wth); // needs re-exec filter->Update() if filter->ReleaseDataFlagOn() https://insightsoftwareconsortium.atlassian.net/browse/ITK-3351?attachmentOrder=desc
    add->SetInput2(th->GetOutput());
    add->ReleaseDataFlagOn();
    add->InPlaceOn(); // overwrites input1
    FilterWatcher watcherA(add);

    typedef itk::SmoothingRecursiveGaussianImageFilter<InputImageType, InputImageType> SmoothType;
    typename SmoothType::Pointer smooth= SmoothType::New();
    smooth->SetInput(add->GetOutput());
    smooth->SetSigma(input->GetSpacing().GetVnlVector().min_value()); // std::min(input->GetSpacing()) not working due to incompatible types  https://github.com/InsightSoftwareConsortium/ITK/blob/master/Modules/ThirdParty/VNL/src/vxl/core/vnl/vnl_vector.h
    smooth->NormalizeAcrossScaleOff();
    smooth->ReleaseDataFlagOn();
    smooth->InPlaceOn();
    FilterWatcher watcherSM(smooth);

    typedef itk::MaskImageFilter<InputImageType, InputImageType, InputImageType> MaskType;
    typename MaskType::Pointer mask = MaskType::New();
    mask->SetInput(smooth->GetOutput());
    mask->SetMaskImage(th->GetOutput()); // needs re-exec th->Update() if th->ReleaseDataFlagOn()
    // mask->ReleaseDataFlagOn(); // needed by stat and mult
    mask->InPlaceOn(); // overwrites SetInput
    FilterWatcher watcherM(mask);

    stat->SetInput(mask->GetOutput());
    stat->Update();

    typedef itk::MultiplyImageFilter<InputImageType, InputImageType, OutputImageType> MultType;
    typename MultType::Pointer mult = MultType::New();
    mult->SetInput(mask->GetOutput());
    mult->SetConstant(1.0 / stat->GetMaximum());
    FilterWatcher watcherMU(mult);

    const typename OutputImageType::Pointer& output= mult->GetOutput();

    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer);
    writer->SetFileName(argv[2]);
    writer->SetInput(mult->GetOutput());
    if(noSDI){
	writer->SetUseCompression(CompChunk);
	}
    else{
	writer->UseCompressionOff(); // writing compressed is not supported when streaming!
	writer->SetNumberOfStreamDivisions(CompChunk);
	}
    try{
        writer->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    return EXIT_SUCCESS;

    }


template<typename InputComponentType, typename InputPixelType>
int dispatch_D(size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    switch (dimensionType){
    case 2:
        res= DoIt<InputComponentType, InputPixelType, 2>(argc, argv);
        break;
    case 3:
        res= DoIt<InputComponentType, InputPixelType, 3>(argc, argv);
        break;
    default:
        std::cerr << "Error: Images of dimension " << dimensionType << " are not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

template<typename InputComponentType>
int dispatch_pT(itk::ImageIOBase::IOPixelType pixelType, size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (pixelType){
    case itk::ImageIOBase::SCALAR:{ // 1 component per pixel
        typedef InputComponentType InputPixelType;
        res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UNKNOWNPIXELTYPE:
    default:
        std::cerr << std::endl << "Error: Pixel type not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

int dispatch_cT(itk::ImageIOBase::IOComponentType componentType, itk::ImageIOBase::IOPixelType pixelType, size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;

    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#a8dc783055a0af6f0a5a26cb080feb178
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00107
    //IOComponentType: UNKNOWNCOMPONENTTYPE, UCHAR, CHAR, USHORT, SHORT, UINT, INT, ULONG, LONG, FLOAT, DOUBLE

    switch (componentType){
    case itk::ImageIOBase::FLOAT:{        // float32
        typedef float InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::DOUBLE:{       // float64
        typedef double InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
    default:
        std::cerr << "unknown component type" << std::endl;
        break;
        }//switch
    return res;
    }


////from http://itk-users.7.n7.nabble.com/Pad-image-with-0-but-keep-its-type-what-ever-it-is-td27442.html
//namespace itk{
  // Description:
  // Get the PixelType and ComponentType from fileName

void GetImageType (std::string fileName,
    itk::ImageIOBase::IOPixelType &pixelType,
    itk::ImageIOBase::IOComponentType &componentType,
    size_t &dimensionType
    ){
    typedef itk::Image<char, 1> ImageType; //template initialization parameters need to be given but can be arbitrary here
    itk::ImageFileReader<ImageType>::Pointer imageReader= itk::ImageFileReader<ImageType>::New();
    imageReader->SetFileName(fileName.c_str());
    imageReader->UpdateOutputInformation();

    if(!imageReader->GetImageIO()->CanStreamRead())
        std::cerr << "Cannot stream the reading of the input. Streaming will be inefficient!" << std::endl;

    pixelType = imageReader->GetImageIO()->GetPixelType();
    componentType = imageReader->GetImageIO()->GetComponentType();
    dimensionType= imageReader->GetImageIO()->GetNumberOfDimensions();

    std::cerr << std::endl << "dimensions: " << dimensionType << std::endl;
    std::cerr << "component type: " << imageReader->GetImageIO()->GetComponentTypeAsString(componentType) << std::endl;
    std::cerr << "component size: " << imageReader->GetImageIO()->GetComponentSize() << std::endl;
    std::cerr << "pixel type (string): " << imageReader->GetImageIO()->GetPixelTypeAsString(imageReader->GetImageIO()->GetPixelType()) << std::endl;
    std::cerr << "pixel type: " << pixelType << std::endl << std::endl;

    }



int main(int argc, char *argv[]){
    if ( argc != 6 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Output_Image"
                  << " compress|stream-chunks"
                  << " bglevel"
                  << " sigma"
                  << std::endl;

        std::cerr << std::endl;
        std::cerr << " no-compress: 0, compress: 1, stream > 1" << std::endl;
        return EXIT_FAILURE;
        }

    int CompChunk= atoi(argv[3]);
    std::cerr << std::endl;
    if(CompChunk == 0){
	std::cerr << "Employing no compression and no streaming." << std::endl;
	}
    else if (CompChunk == 1){
	std::cerr << "Employing compression (streaming not possible then)." << std::endl;
	}
    else if (CompChunk > 1){
	std::cerr << "Employing streaming (compression not possible then)." << std::endl;
	}
    else {
	std::cerr << "compress|stream-chunks must be a positive integer" << std::endl;
        return EXIT_FAILURE;
	}

    itk::ImageIOBase::IOPixelType pixelType;
    typename itk::ImageIOBase::IOComponentType componentType;
    size_t dimensionType;


    try {
        GetImageType(argv[1], pixelType, componentType, dimensionType);
        }//try
    catch( itk::ExceptionObject &excep){
        std::cerr << argv[0] << ": exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
        }

    return dispatch_cT(componentType, pixelType, dimensionType, argc, argv);
    }






