////program for itkROIShiftScaleImageFilter Slice-By-Slice
//01: based on template.cxx


#include <complex>

#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkExtractImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include "filter/self-made/itkROIShiftScaleImageFilter.h"
#include <itkSliceBySliceImageFilter.h>
#include <itkImageFileWriter.h>



template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){

    if( argc != 3 + 2*Dimension){
	fprintf(stderr, "2 + 2*Dimension = %d parameters are needed!\n", 3 + 2*Dimension - 1);
	return EXIT_FAILURE;
	}
	

    typedef InputPixelType  OutputPixelType;

    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;


    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(argv[1]);
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
	
    const typename InputImageType::Pointer& input= reader->GetOutput();


    unsigned int i;
    typename InputImageType::IndexType ROIindex;
    typename InputImageType::SizeType  ROIsize;

    for (i= 0; i < Dimension; i++){
        ROIindex[i]= atoi(argv[3+i]);
    	std::cerr << ROIindex[i] << std::endl;
    	}
    for (i= 0; i < Dimension; i++){
        ROIsize[i]= atoi(argv[3+Dimension+i]);
    	std::cerr << ROIsize[i] << std::endl;
    	}
    typename InputImageType::RegionType ROI(ROIindex, ROIsize);


    typedef itk::ExtractImageFilter<InputImageType, InputImageType> FilterType;
    typename FilterType::Pointer filter= FilterType::New();
    filter->SetInput(input);
    filter->SetExtractionRegion(ROI);
    filter->SetDirectionCollapseToIdentity(); // This is required.

    FilterWatcher watcher1(filter);
    try{
	filter->Update();
	}
    catch(itk::ExceptionObject &ex){
	std::cerr << ex << std::endl;
	return EXIT_FAILURE;
	}

    typedef itk::StatisticsImageFilter<InputImageType> StatType;
    typename StatType::Pointer stat = StatType::New();
    stat->SetInput(filter->GetOutput());
    //stat->GetOutput()->SetRequestedRegion(ROI); //does not work for StatisticsImageFilter
    //stat->GenerateInputRequestedRegion(); //does not help either
    FilterWatcher watcher2(stat);
    try { 
	stat->Update();
        }
    catch(itk::ExceptionObject &ex){ 
	std::cerr << ex << std::endl;
	return EXIT_FAILURE;
	}

    std::cerr << "Min: " << +stat->GetMinimum() << " Max: " << +stat->GetMaximum() << " Mean: " << +stat->GetMean() << " Std: " << +stat->GetSigma() << " Variance: " << +stat->GetVariance() << " Sum: " << +stat->GetSum() << std::endl; //+ promotes variable to a type printable as a number (e.g. for char)

    typename StatType::RealType stackMean= stat->GetMean();
    typename StatType::RealType stackStd= stat->GetSigma();

    typedef itk::SliceBySliceImageFilter<InputImageType, OutputImageType> SBSFilterType;
    typename SBSFilterType::Pointer sbs = SBSFilterType::New();
    sbs->SetInput(input);

    typedef typename  SBSFilterType::InternalInputImageType SBSInputImageType;
    typedef typename  SBSFilterType::InternalOutputImageType SBSOutputImageType;


    typename SBSInputImageType::IndexType ROISBSindex;
    typename SBSInputImageType::SizeType  ROISBSsize;
    for (i= 0; i < Dimension - 1; i++){
        ROISBSindex[i]= atoi(argv[3+i]);
    	std::cerr << ROISBSindex[i] << std::endl;
    	}
    for (i= 0; i < Dimension - 1; i++){
        ROISBSsize[i]= atoi(argv[3+Dimension+i]);
    	std::cerr << ROISBSsize[i] << std::endl;
    	}
    typename SBSInputImageType::RegionType ROISBS(ROISBSindex, ROISBSsize);



    typedef itk::ROIShiftScaleImageFilter<SBSInputImageType, SBSOutputImageType> ROISSType;
    typename ROISSType::Pointer roiss= ROISSType::New();
    //roiss->SetNumberOfThreads(1);//rem. because SBS is NOT multi-threaded
    roiss->SetDesiredMean(stackMean);
    roiss->SetDesiredStd(stackStd);
    roiss->SetROI(ROISBS);

    sbs->SetFilter(roiss);
    FilterWatcher watcher_sbs(sbs);
    try { 
	sbs->Update();
        }
    catch(itk::ExceptionObject &ex){ 
	std::cerr << ex << std::endl;
	return EXIT_FAILURE;
	}

    const typename OutputImageType::Pointer& output= sbs->GetOutput();

    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer);
    writer->SetFileName(argv[2]);
    writer->SetInput(output);
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
    // case itk::ImageIOBase::COMPLEX:{ // 2 components per pixel
    //     typedef std::complex<InputComponentType> InputPixelType;
    //     res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::RGB:{ // 3 components per pixel, limited [0,1]
    //     typedef itk::RGBPixel<InputComponentType> InputPixelType;
    //     res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::RGBA:{ // 4 components per pixel, limited [0,1]
    //     typedef itk::RGBAPixel<InputComponentType> InputPixelType;
    //     res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
    //     } break;
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
    case itk::ImageIOBase::UCHAR:{        // uint8_t
        typedef unsigned char InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::CHAR:{         // int8_t
        typedef char InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::USHORT:{       // uint16_t
        typedef unsigned short InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::SHORT:{        // int16_t
        typedef short InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UINT:{         // uint32_t
        typedef unsigned int InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::INT:{          // int32_t
        typedef int InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::ULONG:{        // uint64_t
        typedef unsigned long InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::LONG:{         // int64_t
        typedef long InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
        } break;
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
    if ( argc < 4 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Output_Image"
		  << " index... size..."
                  << std::endl;

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






