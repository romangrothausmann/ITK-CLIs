////program to use itk to extract a subimage
/// http://www.itk.org/Wiki/ITK/Examples/ImageProcessing/ExtractImageFilter


#include <complex>

#include <itkImageFileReader.h>
#include <itkExtractImageFilter.h>
#include <itkImageFileWriter.h>

#include "itkFilterWatcher.h" 



void dispatch_cT(itk::ImageIOBase::IOPixelType, itk::ImageIOBase::IOComponentType, size_t, int, char **);

template<typename InputComponentType>
void dispatch_pT(itk::ImageIOBase::IOPixelType pixelType, size_t, int, char **);

template<typename InputComponentType, typename InputPixelType>
void dispatch_D(size_t, int, char **);

template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int, char *argv[]);





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
    FilterWatcher watcherI(reader, "reading");
    try{ 
        reader->Update();
        }
    catch(itk::ExceptionObject &ex){ 
	std::cerr << ex << std::endl;
	return EXIT_FAILURE;
	}


    unsigned int i;
    typename InputImageType::IndexType desiredStart;
    typename InputImageType::SizeType desiredSize;

    for (i= 0; i < Dimension; i++){
        desiredStart[i]= atoi(argv[3+i]);
    	std::cerr << desiredStart[i] << std::endl;
    	}
    for (i= 0; i < Dimension; i++){
        desiredSize[i]= atoi(argv[3+Dimension+i]);
    	std::cerr << desiredSize[i] << std::endl;
    	}

    typename InputImageType::RegionType desiredRegion(desiredStart, desiredSize);

    typedef itk::ExtractImageFilter<InputImageType, OutputImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(reader->GetOutput());
    filter->SetExtractionRegion(desiredRegion);
    filter->SetDirectionCollapseToIdentity(); // This is required.

    FilterWatcher watcher1(filter);
    try{ 
        filter->Update();
        }
    catch(itk::ExceptionObject &ex){ 
	std::cerr << ex << std::endl;
	return EXIT_FAILURE;
	}

    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer, "writing");
    writer->SetFileName(argv[2]);
    writer->SetInput(filter->GetOutput());
    //writer->UseCompressionOn();
    try{ 
        writer->Update();
        }
    catch(itk::ExceptionObject &ex){ 
        std::cerr << "itk::ImageFileWriter: " << ex << std::endl;
        return EXIT_FAILURE;
        }

    return EXIT_SUCCESS;

    }


void dispatch_cT(itk::ImageIOBase::IOComponentType componentType, itk::ImageIOBase::IOPixelType pixelType, size_t dimensionType, int argc, char *argv[]){

  //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#a8dc783055a0af6f0a5a26cb080feb178
  //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00107
  //IOComponentType: UNKNOWNCOMPONENTTYPE, UCHAR, CHAR, USHORT, SHORT, UINT, INT, ULONG, LONG, FLOAT, DOUBLE

  switch (componentType){
  case itk::ImageIOBase::UCHAR:{
    typedef unsigned char InputComponentType;
    dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::CHAR:{
    typedef char InputComponentType;
    dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::USHORT:{
    typedef unsigned short InputComponentType;
    dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::SHORT:{
    typedef short InputComponentType;
    dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::UINT:{
    typedef unsigned int InputComponentType;
    dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::INT:{
    typedef int InputComponentType;
    dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::ULONG:{
    typedef unsigned long InputComponentType;
    dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::LONG:{
    typedef long InputComponentType;
    dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::FLOAT:{
    typedef float InputComponentType;
    dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::DOUBLE:{
    typedef double InputComponentType;
    dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
  default:
    std::cout << "unknown component type" << std::endl;
    break;
  }//switch
}

template<typename InputComponentType>
void dispatch_pT(itk::ImageIOBase::IOPixelType pixelType, size_t dimensionType, int argc, char *argv[]){

    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX 

  switch (pixelType){
  case itk::ImageIOBase::SCALAR:{
    typedef InputComponentType InputPixelType;
    dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::RGB:{
    typedef itk::RGBPixel<InputComponentType> InputPixelType;
    dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::RGBA:{
    typedef itk::RGBAPixel<InputComponentType> InputPixelType;
    dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::COMPLEX:{
    typedef std::complex<InputComponentType> InputPixelType;
    dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::VECTOR:{
    typedef itk::VariableLengthVector<InputComponentType> InputPixelType;
    dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::UNKNOWNPIXELTYPE:
  default:
    std::cerr << std::endl << "Error: Pixel type not handled!" << std::endl;
    break;
  }//switch 
}


template<typename InputComponentType, typename InputPixelType>
void dispatch_D(size_t dimensionType, int argc, char *argv[]){
  switch (dimensionType){
  // case 1:
  //   DoIt<InputComponentType, InputPixelType, 1>(argc, argv);
  //   break;
  case 2:
    DoIt<InputComponentType, InputPixelType, 2>(argc, argv);
    break;
  case 3:
    DoIt<InputComponentType, InputPixelType, 3>(argc, argv);
    break;
  default: 
    std::cerr << "Error: Images of dimension " << dimensionType << " are not handled!" << std::endl;
    break;
  }
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
    typedef itk::Image<unsigned char, 3> ImageType;
    itk::ImageFileReader<ImageType>::Pointer imageReader= itk::ImageFileReader<ImageType>::New();
    imageReader->SetFileName(fileName.c_str());
    imageReader->UpdateOutputInformation();

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
	dispatch_cT(componentType, pixelType, dimensionType, argc, argv);
        }//try
    catch( itk::ExceptionObject &excep){
        std::cerr << argv[0] << ": exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
        }
  
    return EXIT_SUCCESS;
    }






