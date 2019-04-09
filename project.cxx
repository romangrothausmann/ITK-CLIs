////program to project an image along an axis (based on itkProjectionImageFilter, which is newer than itkAccumulateImageFilter and itkGetAverageSliceImageFilter?)
//01: based on template_var-filter.cxx


#include <complex>

#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkMeanProjectionImageFilter.h>
#include <itkMedianProjectionImageFilter.h>
#include <itkMaximumProjectionImageFilter.h>
#include <itkMinimumProjectionImageFilter.h>
#include <itkSumProjectionImageFilter.h>
#include <itkStandardDeviationProjectionImageFilter.h>
#include <itkBinaryProjectionImageFilter.h>
#include <itkBinaryThresholdProjectionImageFilter.h>
#include <itkImageFileWriter.h>




template<typename InputComponentType, typename InputPixelType, size_t Dimension, typename InputImageType, typename OutputImageType, typename FilterType>
int DoIt2(int argc, char *argv[], FilterType* filter){

    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(argv[1]);
    reader->ReleaseDataFlagOn();
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

    filter->SetInput(input);
    filter->SetProjectionDimension(atoi(argv[5]));
    filter->ReleaseDataFlagOn();

    FilterWatcher watcher1(filter);
    try{
        filter->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }


    const typename OutputImageType::Pointer& output= filter->GetOutput();

    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer);
    writer->SetFileName(argv[2]);
    writer->SetInput(output);
    writer->SetUseCompression(atoi(argv[3]));
    try{
        writer->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    return EXIT_SUCCESS;

    }


template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){
    int res= 0;

    typedef InputPixelType  OutputPixelType;
    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<OutputPixelType, Dimension - 1>  OutputImageType;


    switch(atoi(argv[4])){
    case 0: {
        typedef itk::MeanProjectionImageFilter<InputImageType, OutputImageType> FilterType;
        typename FilterType::Pointer filter= FilterType::New();
        std::cerr << "Using filter: " << filter->GetNameOfClass() << std::endl;
        res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, OutputImageType, FilterType>(argc, argv, filter);
	} break;
    case 1: {
        typedef itk::MedianProjectionImageFilter<InputImageType, OutputImageType> FilterType;
        typename FilterType::Pointer filter= FilterType::New();
        std::cerr << "Using filter: " << filter->GetNameOfClass() << std::endl;
        res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, OutputImageType, FilterType>(argc, argv, filter);
	} break;
    case 2: {
        typedef itk::MaximumProjectionImageFilter<InputImageType, OutputImageType> FilterType;
        typename FilterType::Pointer filter= FilterType::New();
        std::cerr << "Using filter: " << filter->GetNameOfClass() << std::endl;
        res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, OutputImageType, FilterType>(argc, argv, filter);
	} break;
    case 3: {
        typedef itk::MinimumProjectionImageFilter<InputImageType, OutputImageType> FilterType;
        typename FilterType::Pointer filter= FilterType::New();
        std::cerr << "Using filter: " << filter->GetNameOfClass() << std::endl;
        res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, OutputImageType, FilterType>(argc, argv, filter);
	} break;
    case 4: {
        typedef itk::SumProjectionImageFilter<InputImageType, OutputImageType> FilterType;
        typename FilterType::Pointer filter= FilterType::New();
        std::cerr << "Using filter: " << filter->GetNameOfClass() << std::endl;
        res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, OutputImageType, FilterType>(argc, argv, filter);
	} break;
    case 5: {
        typedef itk::StandardDeviationProjectionImageFilter<InputImageType, OutputImageType> FilterType;
        typename FilterType::Pointer filter= FilterType::New();
        std::cerr << "Using filter: " << filter->GetNameOfClass() << std::endl;
        res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, OutputImageType, FilterType>(argc, argv, filter);
	} break;
    case 6: {
        typedef itk::BinaryProjectionImageFilter<InputImageType, OutputImageType> FilterType;
        typename FilterType::Pointer filter= FilterType::New();
        std::cerr << "Using filter: " << filter->GetNameOfClass() << std::endl;
        res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, OutputImageType, FilterType>(argc, argv, filter);
	} break;
    case 7: {
        typedef itk::BinaryThresholdProjectionImageFilter<InputImageType, OutputImageType> FilterType;
        typename FilterType::Pointer filter= FilterType::New();
        std::cerr << "Using filter: " << filter->GetNameOfClass() << std::endl;
        res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, OutputImageType, FilterType>(argc, argv, filter);
	} break;
    default:
        std::cerr << "unknown filter type." << std::endl;
        res= EXIT_FAILURE;
        break;
        }//switch
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
    case itk::ImageIOBase::SCALAR:{
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
    if ( argc != 6 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Output_Image"
                  << " compress"
                  << " filterType"
                  << " proj-dim"
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






