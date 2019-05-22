////program to use itk to extract a subimage
/// http://www.itk.org/Wiki/ITK/Examples/ImageProcessing/ExtractImageFilter
//02: old version now based on template_vec.cxx
//03: single program for SDI and noSDI (based on new template_vec.cxx)


#include <complex>

#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkExtractImageFilter.h>
#include <itkImageFileWriter.h>



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

    int CompChunk= atoi(argv[3]);
    bool noSDI= CompChunk <= 1; // SDI only if CompChunk > 1

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

    const typename InputImageType::Pointer& input= reader->GetOutput();



    typedef itk::ExtractImageFilter<InputImageType, OutputImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(input);
    filter->SetExtractionRegion(desiredRegion);
    filter->SetDirectionCollapseToIdentity(); // This is required.
    filter->ReleaseDataFlagOn();
    filter->InPlaceOn();

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


    const typename OutputImageType::Pointer& output= filter->GetOutput();

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
    // case 1:
    //     res= dispatch_pT<InputComponentType, CompPerPixel, 1>(pixelType, argc, argv);
    //     break;
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
    if ( argc < 5 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Output_Image"
                  << " compress|stream-chunks"
                  << " index... size..."
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






