////program to use itk to convert between file-formats, streaming version
//01: based on file_converter_SDI.cxx and http://www.itk.org/Wiki/ITK/Examples/IO/VolumeFromSlices
//02: based on template_vec.cxx
//03: single program for SDI and noSDI (based on new template_vec.cxx)


#include <complex>

#include "itkFilterWatcher.h"
#include <itkImageSeriesReader.h>
#include <itkChangeInformationImageFilter.h>
#include <itkImageFileWriter.h>



template<typename InputComponentType, typename InputPixelType, size_t CompPerPixel, size_t Dimension>
int DoIt(int argc, char *argv[]){

    const char offset= 3;
    typedef InputPixelType  OutputPixelType;

    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;

    bool noSDI= false; // SDI only program

    std::vector<std::string> names;

    for(unsigned int i = offset; i < argc; ++i)
        names.push_back(argv[i]);

    // List the files
    for(unsigned int i = 0; i < names.size(); ++i)
        std::cerr << "File: " << names[i] << std::endl;

    typedef itk::ImageSeriesReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    ////reading compressed MHA/MHD is supported for streaming!
    reader->SetFileNames(names);
    reader->ReleaseDataFlagOn();
    if(noSDI){ // not used, SDI only program
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

    std::cerr << "input spacing: " << input->GetSpacing() << std::endl;

    typename InputImageType::SpacingType outputSpacing= input->GetSpacing();
    outputSpacing[Dimension-1]= atof(argv[2]);

    std::cerr << "output spacing: " << outputSpacing << std::endl;


    typedef itk::ChangeInformationImageFilter<InputImageType> FilterType;
    typename FilterType::Pointer filter= FilterType::New();
    filter->SetInput(input);
    filter->ReleaseDataFlagOn();
    filter->SetOutputSpacing(outputSpacing);
    filter->ChangeSpacingOn();

    FilterWatcher watcher1(filter);
    try{
	filter->UpdateOutputInformation();
	}
    catch(itk::ExceptionObject &ex){
	std::cerr << ex << std::endl;
	return EXIT_FAILURE;
	}


    const typename OutputImageType::Pointer& output= filter->GetOutput();

    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer);
    writer->SetFileName(argv[1]);
    writer->SetInput(output);
    writer->UseCompressionOff(); // writing compressed is not supported when streaming!
    writer->SetNumberOfStreamDivisions(names.size());
    try{
        writer->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    return EXIT_SUCCESS;

    }


template<typename InputComponentType, size_t CompPerPixel, size_t Dimension>
int dispatch_pT(itk::ImageIOBase::IOPixelType pixelType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (pixelType){
    case itk::ImageIOBase::SCALAR:{ // 1 component per pixel
        typedef InputComponentType InputPixelType;
        res= DoIt<InputComponentType, InputPixelType, CompPerPixel, Dimension>(argc, argv);
        } break;
    case itk::ImageIOBase::COMPLEX:{ // 2 components per pixel
        typedef std::complex<InputComponentType> InputPixelType;
        res= DoIt<InputComponentType, InputPixelType, CompPerPixel, Dimension>(argc, argv);
        } break;
    case itk::ImageIOBase::RGB:{ // 3 components per pixel, limited [0,1]
        typedef itk::RGBPixel<InputComponentType> InputPixelType;
        res= DoIt<InputComponentType, InputPixelType, CompPerPixel, Dimension>(argc, argv);
        } break;
    case itk::ImageIOBase::RGBA:{ // 4 components per pixel, limited [0,1]
        typedef itk::RGBAPixel<InputComponentType> InputPixelType;
        res= DoIt<InputComponentType, InputPixelType, CompPerPixel, Dimension>(argc, argv);
        } break;
    case itk::ImageIOBase::VECTOR:{
        typedef itk::Vector<InputComponentType, CompPerPixel> InputPixelType;
        res= DoIt<InputComponentType, InputPixelType, CompPerPixel, Dimension>(argc, argv);
        } break;
    case itk::ImageIOBase::UNKNOWNPIXELTYPE:
    default:
        std::cerr << std::endl << "Error: Pixel type not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

template<typename InputComponentType, size_t CompPerPixel>
int dispatch_D(itk::ImageIOBase::IOPixelType pixelType, size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    switch (dimensionType){
    case 1:
        res= dispatch_pT<InputComponentType, CompPerPixel, 1>(pixelType, argc, argv);
        break;
    case 2:
        res= dispatch_pT<InputComponentType, CompPerPixel, 2>(pixelType, argc, argv);
        break;
    case 3:
        res= dispatch_pT<InputComponentType, CompPerPixel, 3>(pixelType, argc, argv);
        break;
    default:
        std::cerr << "Error: Images of dimension " << dimensionType << " are not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

template<typename InputComponentType>
int dispatch_cPP(size_t compPerPixel, itk::ImageIOBase::IOPixelType pixelType, size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    switch (compPerPixel){
    case 1:
        res= dispatch_D<InputComponentType, 1>(pixelType, dimensionType, argc, argv);
        break;
    case 2:
        res= dispatch_D<InputComponentType, 2>(pixelType, dimensionType, argc, argv);
        break;
    case 3:
        res= dispatch_D<InputComponentType, 3>(pixelType, dimensionType, argc, argv);
        break;
    // case 4:
    //     res= dispatch_D<InputComponentType, 4>(pixelType, dimensionType, argc, argv);
    //     break;
    // case 5:
    //     res= dispatch_D<InputComponentType, 5>(pixelType, dimensionType, argc, argv);
    //     break;
    default:
        std::cerr << "Error: NumberOfComponentsPerPixel (" << compPerPixel << ") not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

int dispatch_cT(itk::ImageIOBase::IOComponentType componentType, size_t compPerPixel, itk::ImageIOBase::IOPixelType pixelType, size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;

    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#a8dc783055a0af6f0a5a26cb080feb178
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00107
    //IOComponentType: UNKNOWNCOMPONENTTYPE, UCHAR, CHAR, USHORT, SHORT, UINT, INT, ULONG, LONG, FLOAT, DOUBLE

    switch (componentType){
    case itk::ImageIOBase::UCHAR:{        // uint8_t
        typedef unsigned char InputComponentType;
        res= dispatch_cPP<InputComponentType>(compPerPixel, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::CHAR:{         // int8_t
        typedef char InputComponentType;
        res= dispatch_cPP<InputComponentType>(compPerPixel, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::USHORT:{       // uint16_t
        typedef unsigned short InputComponentType;
        res= dispatch_cPP<InputComponentType>(compPerPixel, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::SHORT:{        // int16_t
        typedef short InputComponentType;
        res= dispatch_cPP<InputComponentType>(compPerPixel, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UINT:{         // uint32_t
        typedef unsigned int InputComponentType;
        res= dispatch_cPP<InputComponentType>(compPerPixel, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::INT:{          // int32_t
        typedef int InputComponentType;
        res= dispatch_cPP<InputComponentType>(compPerPixel, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::ULONG:{        // uint64_t
        typedef unsigned long InputComponentType;
        res= dispatch_cPP<InputComponentType>(compPerPixel, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::LONG:{         // int64_t
        typedef long InputComponentType;
        res= dispatch_cPP<InputComponentType>(compPerPixel, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::FLOAT:{        // float32
        typedef float InputComponentType;
        res= dispatch_cPP<InputComponentType>(compPerPixel, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::DOUBLE:{       // float64
        typedef double InputComponentType;
        res= dispatch_cPP<InputComponentType>(compPerPixel, pixelType, dimensionType, argc, argv);
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
    size_t &compPerPixel,
    size_t &dimensionType
    ){
    typedef itk::VectorImage<char, 1> ImageType; //template initialization parameters need to be given but can be arbitrary here
    itk::ImageFileReader<ImageType>::Pointer imageReader= itk::ImageFileReader<ImageType>::New();
    imageReader->SetFileName(fileName.c_str());
    imageReader->UpdateOutputInformation();

    if(!imageReader->GetImageIO()->CanStreamRead())
        std::cerr << "Cannot stream the reading of the input. Streaming will be inefficient!" << std::endl;

    pixelType = imageReader->GetImageIO()->GetPixelType();
    componentType = imageReader->GetImageIO()->GetComponentType();
    dimensionType= imageReader->GetImageIO()->GetNumberOfDimensions();
    compPerPixel = imageReader->GetOutput()->GetNumberOfComponentsPerPixel(); // needs VectorImage

    std::cerr << std::endl << "dimensions: " << dimensionType << std::endl;
    std::cerr << "component type: " << imageReader->GetImageIO()->GetComponentTypeAsString(componentType) << std::endl;
    std::cerr << "component size: " << imageReader->GetImageIO()->GetComponentSize() << std::endl;
    std::cerr << "pixel type (string): " << imageReader->GetImageIO()->GetPixelTypeAsString(imageReader->GetImageIO()->GetPixelType()) << std::endl;
    std::cerr << "pixel type: " << pixelType << std::endl << std::endl;
    std::cerr << "NumberOfComponentsPerPixel: " << compPerPixel << std::endl;

    }



int main(int argc, char *argv[]){
    if ( argc < 4 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Output_Image"
                  << " spacing-of-last-dim"
                  << " Input_Images"
                  << std::endl;

        return EXIT_FAILURE;
        }

    //// for compression use file_converter afterwards
    std::cerr << "Employing streaming (compression not possible then)." << std::endl;

    itk::ImageIOBase::IOPixelType pixelType;
    typename itk::ImageIOBase::IOComponentType componentType;
    size_t compPerPixel;
    size_t dimensionType;


    try {
        GetImageType(argv[3], pixelType, componentType, compPerPixel, dimensionType);
        }//try
    catch( itk::ExceptionObject &excep){
        std::cerr << argv[0] << ": exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
        }

    return dispatch_cT(componentType, compPerPixel, pixelType, dimensionType+1, argc, argv);
    }






