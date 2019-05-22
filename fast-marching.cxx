////program for itkFastMarchingImageFilter
//01: based on template_2inputs.cxx and snippets from https://code.google.com/p/manageditk/source/browse/trunk/Source/Modules/LevelSetFilters/itkFastMarchingImageFilter.txx?r=2


#include <complex>

#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkFastMarchingImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkImageFileWriter.h>




template<typename InputComponent1, typename TypeInputComponentType2, typename InputPixelType1, typename InputPixelType2, size_t Dimension>
int DoIt(int argc, char *argv[]){

#ifdef USE_FLOAT
    typedef float   OutputPixelType;
    std::cerr << "Using single precision (float)." << std::endl;
#else
    typedef double  OutputPixelType;
    std::cerr << "Using double precision (double)." << std::endl;
#endif

    typedef itk::Image<InputPixelType1, Dimension>  InputImageType1;
    typedef itk::Image<InputPixelType2, Dimension>  InputImageType2;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;


    typedef itk::ImageFileReader<InputImageType1> ReaderType1;
    typename ReaderType1::Pointer reader1 = ReaderType1::New();

    reader1->SetFileName(argv[1]);
    reader1->ReleaseDataFlagOn();
    FilterWatcher watcherI1(reader1);
    watcherI1.QuietOn();
    watcherI1.ReportTimeOn();
    try{
        reader1->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    const typename InputImageType1::Pointer& input1= reader1->GetOutput();


    typedef itk::ImageFileReader<InputImageType2> ReaderType2;
    typename ReaderType2::Pointer reader2 = ReaderType2::New();

    reader2->SetFileName(argv[2]);
    //reader2->ReleaseDataFlagOn();//needed twice
    FilterWatcher watcherI2(reader2);
    watcherI2.QuietOn();
    watcherI2.ReportTimeOn();
    try{
        reader2->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    const typename InputImageType2::Pointer& input2= reader2->GetOutput();

    ////convert input image to speed image, i.e. scale to [0; 1]
    typedef itk::RescaleIntensityImageFilter<InputImageType2, OutputImageType> RescaleType;
    typename RescaleType::Pointer rescaler2 = RescaleType::New();
    rescaler2->SetInput(input2);
    rescaler2->SetOutputMinimum(0);
    rescaler2->SetOutputMaximum(1);
    rescaler2->Update();

    typedef itk::FastMarchingImageFilter<OutputImageType, OutputImageType> FilterType;
    typename FilterType::Pointer filter= FilterType::New();
    filter->SetInput(rescaler2->GetOutput());
    filter->SetStoppingValue(atof(argv[5]));
    filter->ReleaseDataFlagOn();

    typedef typename FilterType::NodeContainer  NodeContainer;
    typedef typename FilterType::NodeType       NodeType;

    typename NodeContainer::Pointer ANodes = NodeContainer::New();
    ANodes->Initialize();

    typedef itk::ImageRegionConstIteratorWithIndex<InputImageType1> IteratorType;
    IteratorType it(input1, input1->GetLargestPossibleRegion() );

    unsigned int count = 0;
    for(it.GoToBegin(); !it.IsAtEnd(); ++it){
        if(it.Get() > 0){
            NodeType node;
            node.SetIndex(it.GetIndex());
            ANodes->InsertElement(count, node);
            count++;
            }
        }
    //filter->SetAlivePoints(ANodes);//alive points are already part of the object, can be omitted?
    filter->SetTrialPoints(ANodes);//trial points are considered for inclusion, eg the layer of pixels around AlivePoints

    FilterWatcher watcher1(filter);
    try{
        filter->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    ////remove all regions with values above FM-StoppingValue
    typedef itk::ThresholdImageFilter<OutputImageType> MFilterType;
    typename MFilterType::Pointer mask = MFilterType::New();
    mask->SetInput(filter->GetOutput());
    mask->ThresholdAbove(filter->GetStoppingValue());
    if (argc > 6)
	mask->SetOutsideValue(atof(argv[6]));
    mask->InPlaceOn();
    FilterWatcher watcherM(mask);

    const typename OutputImageType::Pointer& output= mask->GetOutput();

    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer);
    writer->SetFileName(argv[3]);
    writer->SetInput(output);
    writer->SetUseCompression(atoi(argv[4]));
    try{
        writer->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    return EXIT_SUCCESS;

    }


template<typename InputComponentType1, typename InputComponentType2, typename InputPixelType1, typename InputPixelType2>
int dispatch_D(size_t dimensionType, int argc, char *argv[]){
    int res= 0;
    switch (dimensionType){
    case 2:
        res= DoIt<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2, 2>(argc, argv);
        break;
    case 3:
        res= DoIt<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2, 3>(argc, argv);
        break;
    default:
        std::cerr << "Error: Images of dimension " << dimensionType << " are not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

template<typename InputComponentType1, typename InputComponentType2, typename InputPixelType1>
int dispatch_pT2(itk::ImageIOBase::IOPixelType pixelType2, size_t dimensionType, int argc, char *argv[]){
    int res= 0;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (pixelType2){
    case itk::ImageIOBase::SCALAR:{
        typedef InputComponentType2 InputPixelType2;
        res= dispatch_D<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2>(dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UNKNOWNPIXELTYPE:
    default:
        std::cerr << std::endl << "Error: Pixel type not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

template<typename InputComponentType1, typename InputComponentType2>
int dispatch_pT1(itk::ImageIOBase::IOPixelType pixelType1, itk::ImageIOBase::IOPixelType pixelType2, size_t dimensionType, int argc, char *argv[]){
    int res= 0;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (pixelType1){
    case itk::ImageIOBase::SCALAR:{
        typedef InputComponentType1 InputPixelType1;
        res= dispatch_pT2<InputComponentType1, InputComponentType2, InputPixelType1>(pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UNKNOWNPIXELTYPE:
    default:
        std::cerr << std::endl << "Error: Pixel type not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

template<typename InputComponentType1>
int dispatch_cT2(itk::ImageIOBase::IOComponentType componentType2, itk::ImageIOBase::IOPixelType pixelType1, itk::ImageIOBase::IOPixelType pixelType2, size_t dimensionType, int argc, char *argv[]){
    int res= 0;

    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#a8dc783055a0af6f0a5a26cb080feb178
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00107
    //IOComponentType: UNKNOWNCOMPONENTTYPE, UCHAR, CHAR, USHORT, SHORT, UINT, INT, ULONG, LONG, FLOAT, DOUBLE

    switch (componentType2){
    case itk::ImageIOBase::UCHAR:{        // uint8_t
        typedef unsigned char InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::CHAR:{         // int8_t
        typedef char InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::USHORT:{       // uint16_t
        typedef unsigned short InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::SHORT:{        // int16_t
        typedef short InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UINT:{         // uint32_t
        typedef unsigned int InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::INT:{          // int32_t
        typedef int InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::ULONG:{        // uint64_t
        typedef unsigned long InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::LONG:{         // int64_t
        typedef long InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::FLOAT:{        // float32
        typedef float InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::DOUBLE:{       // float64
        typedef double InputComponentType2;
        res= dispatch_pT1<InputComponentType1, InputComponentType2>(pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
    default:
        std::cerr << "unknown component type" << std::endl;
        break;
        }//switch
    return res;
    }

int dispatch_cT1(itk::ImageIOBase::IOComponentType componentType1, itk::ImageIOBase::IOComponentType componentType2, itk::ImageIOBase::IOPixelType pixelType1, itk::ImageIOBase::IOPixelType pixelType2, size_t dimensionType, int argc, char *argv[]){
    int res= 0;

    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#a8dc783055a0af6f0a5a26cb080feb178
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00107
    //IOComponentType: UNKNOWNCOMPONENTTYPE, UCHAR, CHAR, USHORT, SHORT, UINT, INT, ULONG, LONG, FLOAT, DOUBLE

    switch (componentType1){
    case itk::ImageIOBase::UCHAR:{        // uint8_t
        typedef unsigned char InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
/*
    case itk::ImageIOBase::CHAR:{         // int8_t
        typedef char InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::USHORT:{       // uint16_t
        typedef unsigned short InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::SHORT:{        // int16_t
        typedef short InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UINT:{         // uint32_t
        typedef unsigned int InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::INT:{          // int32_t
        typedef int InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::ULONG:{        // uint64_t
        typedef unsigned long InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::LONG:{         // int64_t
        typedef long InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::FLOAT:{        // float32
        typedef float InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::DOUBLE:{       // float64
        typedef double InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
*/
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
    if ( argc < 6 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Source_Image"
                  << " Speed_Image"
                  << " Output_Image"
                  << " compress"
                  << " stop-value"
                  << " [not-reached-value]"
                  << std::endl;

        return EXIT_FAILURE;
        }

    itk::ImageIOBase::IOPixelType pixelType1;
    itk::ImageIOBase::IOPixelType pixelType2;
    typename itk::ImageIOBase::IOComponentType componentType1;
    typename itk::ImageIOBase::IOComponentType componentType2;
    size_t dimensionType1;
    size_t dimensionType2;


    try {
        GetImageType(argv[1], pixelType1, componentType1, dimensionType1);
        GetImageType(argv[2], pixelType2, componentType2, dimensionType2);
        }//try
    catch( itk::ExceptionObject &excep){
        std::cerr << argv[0] << ": exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
        }

    if (dimensionType1 != dimensionType2){
        std::cout << "Input images need to be of the same dimension." << std::endl;
        return EXIT_FAILURE;
        }

    return dispatch_cT1(componentType1, componentType2, pixelType1, pixelType2, dimensionType1, argc, argv);
    }






