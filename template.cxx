////program for
//01: based on template.cxx


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



template<typename InputComponentType, size_t NumComponents, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){

    typedef   OutputPixelType;

    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;

    // itk::CStyleCommand::Pointer eventCallbackITK;
    // eventCallbackITK = itk::CStyleCommand::New();
    // eventCallbackITK->SetCallback(FilterEventHandlerITK<InputImageType, OutputImageType>);

    int CompChunk= atoi(argv[3]);
    bool noSDI= CompChunk <= 1; // SDI only if CompChunk > 1

    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(argv[1]);
    reader->ReleaseDataFlagOn();
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



    typedef itk::<InputImageType> FilterType;
    typename FilterType::Pointer filter= FilterType::New();
    filter->SetInput(input);
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
    writer->SetFileName(argv[2]);
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


template<typename InputComponentType, size_t NumComponents, typename InputPixelType>
int dispatch_D(size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    switch (dimensionType){
    case 1:
        res= DoIt<InputComponentType, NumComponents, InputPixelType, 1>(argc, argv);
        break;
    case 2:
        res= DoIt<InputComponentType, NumComponents, InputPixelType, 2>(argc, argv);
        break;
    case 3:
        res= DoIt<InputComponentType, NumComponents, InputPixelType, 3>(argc, argv);
        break;
    default:
        std::cerr << "Error: Images of dimension " << dimensionType << " are not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

template<typename InputComponentType>
int dispatch_pT(size_t numComponentsType, itk::ImageIOBase::IOPixelType pixelType, size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (numComponentsType){ // ignoring pixelType (might still be needed in case of dep. on limited value range [0,1])
    case 1: { // e.g. scalar
	typedef itk::Vector<InputComponentType, 1> InputPixelType;
	res= dispatch_D<InputComponentType, 1, InputPixelType>(dimensionType, argc, argv);
	} break;
    case 2: { // e.g. complex
	typedef itk::Vector<InputComponentType, 2> InputPixelType;
	res= dispatch_D<InputComponentType, 2, InputPixelType>(dimensionType, argc, argv);
	} break;
    case 3: { // e.g. RGB
	typedef itk::Vector<InputComponentType, 3> InputPixelType;
	res= dispatch_D<InputComponentType, 3, InputPixelType>(dimensionType, argc, argv);
	} break;
    case 4: { // e.g. RGBA
	typedef itk::Vector<InputComponentType, 4> InputPixelType;
	res= dispatch_D<InputComponentType, 4, InputPixelType>(dimensionType, argc, argv);
	} break;
    default:
	std::cerr << std::endl << "Error: Number of components (Pixel type) not handled!" << std::endl;
	break;
	} // numComponentsType switch
    return res;
    }

int dispatch_cT(itk::ImageIOBase::IOComponentType componentType, size_t numComponentsType, itk::ImageIOBase::IOPixelType pixelType, size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;

    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#a8dc783055a0af6f0a5a26cb080feb178
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00107
    //IOComponentType: UNKNOWNCOMPONENTTYPE, UCHAR, CHAR, USHORT, SHORT, UINT, INT, ULONG, LONG, FLOAT, DOUBLE

    switch (componentType){
    case itk::ImageIOBase::UCHAR:{        // uint8_t
        typedef unsigned char InputComponentType;
        res= dispatch_pT<InputComponentType>(numComponentsType, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::CHAR:{         // int8_t
        typedef char InputComponentType;
        res= dispatch_pT<InputComponentType>(numComponentsType, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::USHORT:{       // uint16_t
        typedef unsigned short InputComponentType;
        res= dispatch_pT<InputComponentType>(numComponentsType, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::SHORT:{        // int16_t
        typedef short InputComponentType;
        res= dispatch_pT<InputComponentType>(numComponentsType, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::UINT:{         // uint32_t
        typedef unsigned int InputComponentType;
        res= dispatch_pT<InputComponentType>(numComponentsType, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::INT:{          // int32_t
        typedef int InputComponentType;
        res= dispatch_pT<InputComponentType>(numComponentsType, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::ULONG:{        // uint64_t
        typedef unsigned long InputComponentType;
        res= dispatch_pT<InputComponentType>(numComponentsType, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::LONG:{         // int64_t
        typedef long InputComponentType;
        res= dispatch_pT<InputComponentType>(numComponentsType, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::FLOAT:{        // float32
        typedef float InputComponentType;
        res= dispatch_pT<InputComponentType>(numComponentsType, pixelType, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::DOUBLE:{       // float64
        typedef double InputComponentType;
        res= dispatch_pT<InputComponentType>(numComponentsType, pixelType, dimensionType, argc, argv);
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
    size_t &numComponentsType, 
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
    numComponentsType = imageReader->GetImageIO()->GetNumberOfComponents();
    dimensionType= imageReader->GetImageIO()->GetNumberOfDimensions();

    std::cerr << std::endl << "dimensions: " << dimensionType << std::endl;
    std::cerr << "component type: " << imageReader->GetImageIO()->GetComponentTypeAsString(componentType) << std::endl;
    std::cerr << "component size: " << imageReader->GetImageIO()->GetComponentSize() << std::endl;
    std::cerr << "component #: " << numComponentsType << std::endl;
    std::cerr << "pixel type (string): " << imageReader->GetImageIO()->GetPixelTypeAsString(imageReader->GetImageIO()->GetPixelType()) << std::endl;
    std::cerr << "pixel type: " << pixelType << std::endl << std::endl;

    }



int main(int argc, char *argv[]){
    if ( argc != 4 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Output_Image"
                  << " compress|stream-chunks"
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
    size_t numComponentsType;
    size_t dimensionType;


    try {
        GetImageType(argv[1], pixelType, componentType, numComponentsType, dimensionType);
        }//try
    catch( itk::ExceptionObject &excep){
        std::cerr << argv[0] << ": exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
        }

    return dispatch_cT(componentType, numComponentsType, pixelType, dimensionType, argc, argv);
    }






