////program to mask an image within a ROI, i.e. past the ROI into a black image (different to itkRegionOfInterestImageFilter)
//01: based on subimage_extract.cxx


#include <complex>

#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkPasteImageFilter.h>
#include <itkImageFileWriter.h>

#ifdef USE_SDI
#include <itkPipelineMonitorImageFilter.h>
#endif



template<typename InputComponentType, typename InputPixelType, size_t CompPerPixel, size_t Dimension>
int DoIt(int argc, char *argv[]){

    const char offset= 4;
    if( argc != offset + 2*Dimension){
        fprintf(stderr, "2 + 2*Dimension = %d parameters are needed!\n", offset + 2*Dimension - 1);
        return EXIT_FAILURE;
        }


    typedef InputPixelType  OutputPixelType;

    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;


    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(argv[1]);
    reader->ReleaseDataFlagOn();
    reader->UpdateOutputInformation();

    std::cerr << "input region index: " << reader->GetOutput()->GetLargestPossibleRegion().GetIndex()
              << "  size: " <<  reader->GetOutput()->GetLargestPossibleRegion().GetSize()
              << std::endl;

    unsigned int i;
    typename InputImageType::IndexType desiredStart;
    typename InputImageType::SizeType desiredSize;

    for (i= 0; i < Dimension; i++)
        desiredStart[i]= atoi(argv[offset+i]);
    for (i= 0; i < Dimension; i++)
        desiredSize[i]=  atoi(argv[offset+Dimension+i]);

    typename InputImageType::RegionType desiredRegion(desiredStart, desiredSize);
    std::cerr << "desired region index: " << desiredRegion.GetIndex()
              << "  size: " << desiredRegion.GetSize()
              << std::endl;

    if(!reader->GetOutput()->GetLargestPossibleRegion().IsInside(desiredRegion)){
        std::cerr << "desired region is not inside the largest possible input region! Forgot index -1?" << std::endl;
        return EXIT_FAILURE;
        }
    reader->GetOutput()->SetRequestedRegion(desiredRegion);

    typename InputImageType::Pointer black= InputImageType::New();
    black->SetRegions(reader->GetOutput()->GetLargestPossibleRegion());
    black->Allocate();
    black->FillBuffer(itk::NumericTraits<InputPixelType>::Zero);

#ifndef USE_SDI
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
#endif


    typedef itk::PasteImageFilter<InputImageType, InputImageType, OutputImageType> FilterType;
    typename FilterType::Pointer filter= FilterType::New();
    filter->SetSourceImage(reader->GetOutput());
    filter->SetDestinationImage(black);
    filter->ReleaseDataFlagOn();
    filter->InPlaceOn();
    filter->SetSourceRegion(desiredRegion);
    filter->SetDestinationIndex(desiredRegion.GetIndex());

#ifndef USE_SDI
    FilterWatcher watcher1(filter);
    try{
        filter->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

#else
    typedef itk::PipelineMonitorImageFilter<OutputImageType> MonitorFilterType;
    typename MonitorFilterType::Pointer monitorFilter = MonitorFilterType::New();
    monitorFilter->SetInput(filter->GetOutput());
#endif

    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer);
    writer->SetFileName(argv[2]);
#ifndef USE_SDI
    writer->SetInput(filter->GetOutput());
    writer->SetUseCompression(atoi(argv[3]));
#else
    writer->SetInput(monitorFilter->GetOutput());
    writer->UseCompressionOff(); //writing compressed is not supported for streaming!
    writer->SetNumberOfStreamDivisions(atoi(argv[3]));
#endif
    try{
        writer->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

#ifdef USE_SDI
    if (!monitorFilter->VerifyAllInputCanStream(atoi(argv[3]))){ // reports a warning if expected and actual # chunks differ
        // std::cerr << monitorFilter;
        }
#endif

    return EXIT_SUCCESS;

    }


template<typename InputComponentType, typename InputPixelType, size_t CompPerPixel>
int dispatch_D(size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    switch (dimensionType){
    // case 1:
    //     res= DoIt<InputComponentType, InputPixelType, CompPerPixel, 1>(argc, argv);
    //     break;
    case 2:
        res= DoIt<InputComponentType, InputPixelType, CompPerPixel, 2>(argc, argv);
        break;
    case 3:
        res= DoIt<InputComponentType, InputPixelType, CompPerPixel, 3>(argc, argv);
        break;
    default:
        std::cerr << "Error: Images of dimension " << dimensionType << " are not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

template<typename InputComponentType, size_t CompPerPixel>
int dispatch_pT(itk::ImageIOBase::IOPixelType pixelType, size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (pixelType){
    case itk::ImageIOBase::SCALAR:{
        typedef InputComponentType InputPixelType;
        res= dispatch_D<InputComponentType, InputPixelType, CompPerPixel>(dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::RGB:{
        typedef itk::RGBPixel<InputComponentType> InputPixelType;
        res= dispatch_D<InputComponentType, InputPixelType, CompPerPixel>(dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::RGBA:{
        typedef itk::RGBAPixel<InputComponentType> InputPixelType;
        res= dispatch_D<InputComponentType, InputPixelType, CompPerPixel>(dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::COMPLEX:{
        typedef std::complex<InputComponentType> InputPixelType;
        res= dispatch_D<InputComponentType, InputPixelType, CompPerPixel>(dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::VECTOR:{
        typedef itk::Vector<InputComponentType, CompPerPixel> InputPixelType;
        res= dispatch_D<InputComponentType, InputPixelType, CompPerPixel>(dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UNKNOWNPIXELTYPE:
    default:
        std::cerr << std::endl << "Error: Pixel type not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

template<typename InputComponentType>
int dispatch_cPP(size_t compPerPixel, itk::ImageIOBase::IOPixelType pixelType, size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    switch (compPerPixel){
    case 1:
        res= dispatch_pT<InputComponentType, 1>(pixelType, dimensionType, argc, argv);
        break;
    case 2:
        res= dispatch_pT<InputComponentType, 2>(pixelType, dimensionType, argc, argv);
        break;
    case 3:
        res= dispatch_pT<InputComponentType, 3>(pixelType, dimensionType, argc, argv);
        break;
    // case 4:
    //     res= dispatch_pT<InputComponentType, 4>(pixelType, dimensionType, argc, argv);
    //     break;
    // case 5:
    //     res= dispatch_pT<InputComponentType, 5>(pixelType, dimensionType, argc, argv);
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
    if ( argc < 5 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Output_Image"
#ifndef USE_SDI
                  << " compress"
#else
                  << " stream-chunks"
#endif
                  << " index... size..."
                  << std::endl;

        return EXIT_FAILURE;
        }

    itk::ImageIOBase::IOPixelType pixelType;
    typename itk::ImageIOBase::IOComponentType componentType;
    size_t compPerPixel;
    size_t dimensionType;


    try {
        GetImageType(argv[1], pixelType, componentType, compPerPixel, dimensionType);
        }//try
    catch( itk::ExceptionObject &excep){
        std::cerr << argv[0] << ": exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
        }

    return dispatch_cT(componentType, compPerPixel, pixelType, dimensionType, argc, argv);
    }






