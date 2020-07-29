////program for
//01: based on template_2inputs.cxx


#include <complex>

#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>



// template<typename ReaderImageType, typename WriterImageType>
// void FilterEventHandlerITK(itk::Object *caller, const itk::EventObject &event, void*){

//     const itk::ProcessObject* filter = static_cast<const itk::ProcessObject*>(caller);

//     if(itk::ProgressEvent().CheckEvent(&event))
//         fprintf(stderr, "\r%s progress: %5.1f%%", filter->GetNameOfClass(), 100.0 * filter->GetProgress());//stderr is flushed directly
//     else if(itk::StartEvent().CheckEvent(&event)){
// 	if(strstr(filter->GetNameOfClass(), "ImageFileReader"))
// 	    std::cerr << "Reading: " << (dynamic_cast<itk::ImageFileReader<ReaderImageType> *>(caller))->GetFileName() << std::endl;//cast only works if reader was instanciated for ReaderImageType!
// 	else if(strstr(filter->GetNameOfClass(), "ImageFileWriter"))
// 	    std::cerr << "Writing: " << (dynamic_cast<itk::ImageFileWriter<WriterImageType> *>(caller))->GetFileName() << std::endl;//cast only works if writer was instanciated for WriterImageType!
// 	}
//     else if(itk::IterationEvent().CheckEvent(&event))
//         std::cerr << " Iteration: " << (dynamic_cast<itk::SliceBySliceImageFilter<ReaderImageType, WriterImageType> *>(caller))->GetSliceIndex() << std::endl;
//     else if(itk::EndEvent().CheckEvent(&event))
//         std::cerr << std::endl;
//     }



template<typename InputComponentType1, typename InputComponentType2, typename InputPixelType1, typename InputPixelType2, size_t Dimension>
int DoIt(int argc, char *argv[]){

    typedef   OutputPixelType;

    typedef itk::Image<InputPixelType1, Dimension>  InputImageType1;
    typedef itk::Image<InputPixelType2, Dimension>  InputImageType2;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;

    // itk::CStyleCommand::Pointer eventCallbackITK;
    // eventCallbackITK = itk::CStyleCommand::New();
    // eventCallbackITK->SetCallback(FilterEventHandlerITK<InputImageType, OutputImageType>);


    int CompChunk= atoi(argv[4]);
    bool noSDI= CompChunk <= 1; // SDI only if CompChunk > 1

    typedef itk::ImageFileReader<InputImageType1> ReaderType1;
    typename ReaderType1::Pointer reader1 = ReaderType1::New();

    reader1->SetFileName(argv[1]);
    reader1->ReleaseDataFlagOn();
    if(noSDI){
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
	}
    else{
	reader1->UpdateOutputInformation();
	}

    const typename InputImageType1::Pointer& input1= reader1->GetOutput();


    typedef itk::ImageFileReader<InputImageType2> ReaderType2;
    typename ReaderType2::Pointer reader2 = ReaderType2::New();

    reader2->SetFileName(argv[2]);
    reader2->ReleaseDataFlagOn();
    if(noSDI){
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
	}
    else{
	reader2->UpdateOutputInformation();
	}

    const typename InputImageType2::Pointer& input2= reader2->GetOutput();


    typedef itk::<InputImageType1, InputImageType2, OutputImageType> FilterType;
    typename FilterType::Pointer filter= FilterType::New();
    filter->SetInput1(reader1->GetOutput());
    filter->SetInput2(reader2->GetOutput());
    filter->ReleaseDataFlagOn();
    filter->InPlaceOn();

    if(noSDI){
	FilterWatcher watcher1(filter);
	// filter->AddObserver(itk::ProgressEvent(), eventCallbackITK);
	// filter->AddObserver(itk::IterationEvent(), eventCallbackITK);
	// filter->AddObserver(itk::EndEvent(), eventCallbackITK);
	try{
	    filter->Update();
	    }
	catch(itk::ExceptionObject &ex){
	    std::cerr << ex << std::endl;
	    return EXIT_FAILURE;
	    }
	}


    const typename OutputImageType::Pointer& output= filterXYZ->GetOutput();

    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer);
    writer->SetFileName(argv[3]);
    writer->SetInput(output);
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


template<typename InputComponentType1, typename InputComponentType2, typename InputPixelType1, typename InputPixelType2>
int dispatch_D(size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    switch (dimensionType){
    case 1:
        res= DoIt<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2, 1>(argc, argv);
        break;
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
    int res= EXIT_FAILURE;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (pixelType2){
    case itk::ImageIOBase::SCALAR:{ // 1 component per pixel
        typedef InputComponentType2 InputPixelType2;
        res= dispatch_D<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2>(dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::COMPLEX:{ // 2 components per pixel
        typedef std::complex<InputComponentType2> InputPixelType2;
        res= dispatch_D<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2>(dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::RGB:{ // 3 components per pixel, limited [0,1]
        typedef itk::RGBPixel<InputComponentType2> InputPixelType2;
        res= dispatch_D<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2>(dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::RGBA:{ // 4 components per pixel, limited [0,1]
        typedef itk::RGBAPixel<InputComponentType2> InputPixelType2;
        res= dispatch_D<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2>(dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::VECTOR:{
        typedef itk::VariableLengthVector<InputComponentType2> InputPixelType2;
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
    int res= EXIT_FAILURE;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (pixelType1){
    case itk::ImageIOBase::SCALAR:{ // 1 component per pixel
        typedef InputComponentType1 InputPixelType1;
        res= dispatch_pT2<InputComponentType1, InputComponentType2, InputPixelType1>(pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::COMPLEX:{ // 2 components per pixel
        typedef std::complex<InputComponentType1> InputPixelType1;
        res= dispatch_pT2<InputComponentType1, InputComponentType2, InputPixelType1>(pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::RGB:{ // 3 components per pixel, limited [0,1]
        typedef itk::RGBPixel<InputComponentType1> InputPixelType1;
        res= dispatch_pT2<InputComponentType1, InputComponentType2, InputPixelType1>(pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::RGBA:{ // 4 components per pixel, limited [0,1]
        typedef itk::RGBAPixel<InputComponentType1> InputPixelType1;
        res= dispatch_pT2<InputComponentType1, InputComponentType2, InputPixelType1>(pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::VECTOR:{
        typedef itk::VariableLengthVector<InputComponentType1> InputPixelType1;
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
    int res= EXIT_FAILURE;

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
    int res= EXIT_FAILURE;

    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#a8dc783055a0af6f0a5a26cb080feb178
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00107
    //IOComponentType: UNKNOWNCOMPONENTTYPE, UCHAR, CHAR, USHORT, SHORT, UINT, INT, ULONG, LONG, FLOAT, DOUBLE

    switch (componentType1){
    case itk::ImageIOBase::UCHAR:{        // uint8_t
        typedef unsigned char InputComponentType1;
        res= dispatch_cT2<InputComponentType1>(componentType2, pixelType1, pixelType2, dimensionType, argc, argv);
        } break;
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
                  << " Input_Image1"
                  << " Input_Image2"
                  << " Output_Image"
                  << " compress|stream-chunks"
                  << std::endl;

        std::cerr << std::endl;
        std::cerr << " no-compress: 0, compress: 1, stream > 1" << std::endl;
        return EXIT_FAILURE;
        }

    int CompChunk= atoi(argv[4]);
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






