////program for itkBinaryFillholeImageFilter
//01: based on template.cxx
//02: itkBinaryFillholeImageFilter is based on ShapeOpeningLabelMapFilter and LabelMapMaskImageFilter (and correctly handles images containing more holes than the image-type can hold)
//    basically took contents of BinaryFillholeImageFilter to allow reporting the # of filled holes

#include <iomanip>

#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkBinaryNotImageFilter.h>
#include <itkBinaryImageToShapeLabelMapFilter.h>
#include <itkShapeOpeningLabelMapFilter.h>//takes a LabelMap as input wheras LabelShapeOpeningImageFilter takes a labeled image as input, LabelImageToLabelMapFilter converts such an image into a LabelMap
#include <itkLabelMapMaskImageFilter.h>
#include <itkImageFileWriter.h>




template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){

    typedef InputPixelType  OutputPixelType;

    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;


    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(argv[1]);
    //reader->ReleaseDataFlagOn();//input needed twice: notInput->SetInput(input) and binarizer->SetFeatureImage(input)
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

    ////taken from Modules/Filtering/LabelMap/include/itkBinaryFillholeImageFilter.hxx
    int m_FullyConnected= atoi(argv[4]);
    InputPixelType m_ForegroundValue = itk::NumericTraits<InputPixelType>::max();
    InputPixelType backgroundValue = itk::NumericTraits<InputPixelType>::ZeroValue();
    if ( m_ForegroundValue == backgroundValue )
        {
        // current background value is already used for foreground value
        // choose another one
        backgroundValue = itk::NumericTraits<InputPixelType>::max();
        }

    typedef itk::BinaryNotImageFilter< InputImageType > NotType;
    typename NotType::Pointer notInput = NotType::New();
    notInput->SetInput(input);
    notInput->SetForegroundValue( m_ForegroundValue );
    notInput->SetBackgroundValue( backgroundValue );
    notInput->SetReleaseDataFlag( true );
    FilterWatcher watcher0(notInput);

    typedef typename itk::BinaryImageToShapeLabelMapFilter< InputImageType > LabelizerType;
    typename LabelizerType::Pointer labelizer = LabelizerType::New();
    labelizer->SetInput( notInput->GetOutput() );
    labelizer->SetInputForegroundValue( m_ForegroundValue );
    labelizer->SetOutputBackgroundValue( backgroundValue );
    labelizer->SetFullyConnected( m_FullyConnected );
    FilterWatcher watcher1(labelizer);

    typedef typename LabelizerType::OutputImageType                  LabelMapType;
    typedef typename itk::ShapeOpeningLabelMapFilter< LabelMapType > OpeningType;
    typename OpeningType::Pointer opening= OpeningType::New();
    opening->SetInput(labelizer->GetOutput());
    opening->SetAttribute(LabelMapType::LabelObjectType::NUMBER_OF_PIXELS_ON_BORDER);
    opening->SetLambda(1);
    FilterWatcher watcher2(opening);

    // invert the image during the binarization
    typedef typename itk::LabelMapMaskImageFilter< LabelMapType, OutputImageType > BinarizerType;
    typename BinarizerType::Pointer binarizer = BinarizerType::New();
    binarizer->SetInput(opening->GetOutput());
    binarizer->SetLabel(backgroundValue);
    binarizer->SetNegated(true);
    binarizer->SetBackgroundValue(m_ForegroundValue);
    binarizer->SetFeatureImage(input);//will cause another read if reader->ReleaseDataFlagOn();
    FilterWatcher watcher3(binarizer);
    try{
        binarizer->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    const typename OutputImageType::Pointer& output= binarizer->GetOutput();

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

    //// ShapeOpeningLabelMapFilter calls std::map::erase(key) which seems not only to delete the contents but also reduces the total count of elements, so GetNumberOfLabelObjects() reflects the actually contained elements even if not labelled consecutively
    std::cerr << "# of holes removed: " << labelizer->GetOutput()->GetNumberOfLabelObjects() - opening->GetOutput()->GetNumberOfLabelObjects() << std::endl;

    //// report sizes of filled holes
    opening->SetLambda(0);
    opening->ReverseOrderingOn();
    opening->Update();//has to be placed after writer->Update()!
    for(unsigned int i = 0; i < opening->GetOutput()->GetNumberOfLabelObjects(); ++i){
	std::cerr << "hole: " << std::setw(6) << +i << " size: " << +opening->GetOutput()->GetNthLabelObject(i)->GetNumberOfPixels() << std::endl; //::SizeValueType
	}

    return EXIT_SUCCESS;

    }


template<typename InputComponentType, typename InputPixelType>
int dispatch_D(size_t dimensionType, int argc, char *argv[]){
    int res= 0;
    switch (dimensionType){
    // case 1:
    //     res= DoIt<InputComponentType, InputPixelType, 1>(argc, argv);
    //     break;
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
    int res= 0;
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
    int res= 0;

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
    if ( argc != 5 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Output_Image"
                  << " compress"
                  << " fully-connected"
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






