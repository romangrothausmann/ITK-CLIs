/////program to watershed an image
//02: multi dim, multi type

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkMorphologicalWatershedImageFilter.h>
#include "itkFilterWatcher.h"



template<typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){

    //typedef InputPixelType  OutputPixelType;
    typedef  unsigned short OutputPixelType;
    
    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;

    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
 
    reader->SetFileName(argv[1]);
    FilterWatcher watcherI(reader, "reading");
    try
        { 
        reader->Update();
        }
    catch (itk::ExceptionObject &ex)
        { 
	if (!strcmp(ex.GetDescription(), "Filter does not have progress.")){
	  std::cerr << ex << std::endl;
	  return EXIT_FAILURE;
	  }
        }

    typedef itk::MorphologicalWatershedImageFilter<InputImageType, OutputImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(reader->GetOutput());

    filter->SetLevel(atoi(argv[3]));
    filter->SetFullyConnected(atoi(argv[4]));
    filter->SetMarkWatershedLine(atoi(argv[5]));

    FilterWatcher watcher(filter);
    try { 
        filter->Update();
        }
    catch (itk::ExceptionObject &ex)
        { 
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
        }


    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer, "writing");
    writer->SetFileName(argv[2]);
    writer->SetInput(filter->GetOutput());
    writer->UseCompressionOn();
    try
        { 
        writer->Update();
        }
    catch (itk::ExceptionObject &ex)
        { 
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
        }

    return EXIT_SUCCESS;

    }


////from http://itk-users.7.n7.nabble.com/Pad-image-with-0-but-keep-its-type-what-ever-it-is-td27442.html
//namespace itk{
  // Description:
  // Get the PixelType and ComponentType from fileName

void GetImageType (std::string fileName,
    itk::ImageIOBase::IOPixelType &pixelType,
    itk::ImageIOBase::IOComponentType &componentType,
    size_t &dimensionType
    //ImageIOBase::IODimensionType &dimensionType
    ){
    typedef itk::Image<unsigned char, 3> ImageType;
    itk::ImageFileReader<ImageType>::Pointer imageReader= itk::ImageFileReader<ImageType>::New();
    imageReader->SetFileName(fileName.c_str());
    imageReader->UpdateOutputInformation();

    pixelType = imageReader->GetImageIO()->GetPixelType();
    componentType = imageReader->GetImageIO()->GetComponentType();
    dimensionType= imageReader->GetImageIO()->GetNumberOfDimensions();

    //std::cout << "Pixel Type is " << imageReader->GetImageIO()->GetComponentTypeAsString(pixelType) << std::endl;
    std::cout << "Pixel Type is " << imageReader->GetImageIO()->GetComponentTypeAsString(componentType) << std::endl;
    std::cout << "numDimensions: " << dimensionType << std::endl;  
    std::cout << "component size: " << imageReader->GetImageIO()->GetComponentSize() << std::endl; 
    std::cout << "pixel type (string): " << imageReader->GetImageIO()->GetPixelTypeAsString(imageReader->GetImageIO()->GetPixelType()) << std::endl; 
    std::cout << "pixel type: " << pixelType << std::endl; 

    }
  


int main(int argc, char *argv[]){
    if ( argc != 6 )
	{
	std::cerr << "Missing Parameters: "
		  << argv[0]
		  << " Input_Image"
		  << " Output_Image"
		  << " level connectivity lines"
    		  << std::endl;

	return EXIT_FAILURE;
	}

    std::string ifn = argv[1]; 
    // std::string ofn = argv[2]; 
    // int compress= atoi(argv[3]);

    itk::ImageIOBase::IOPixelType pixelType;
    itk::ImageIOBase::IOComponentType componentType;
    //itk::ImageIOBase::IODimensionType dimensionType;
    size_t dimensionType;

  try
    {
      GetImageType(ifn, pixelType, componentType, dimensionType);

      ////switch dim first as dim=1 does not work for MorphologicalWatershedImageFilter

      switch (dimensionType){
      case 2:
	switch (componentType){
	case itk::ImageIOBase::UCHAR:{
	  typedef unsigned char InputPixelType;
	  DoIt<InputPixelType, 2>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::CHAR:{
	  typedef char InputPixelType;
	  DoIt<InputPixelType, 2>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::USHORT:{
	  typedef unsigned short InputPixelType;
	  DoIt<InputPixelType, 2>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::SHORT:{
	  typedef short InputPixelType;
	  DoIt<InputPixelType, 2>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::UINT:{
	  typedef unsigned int InputPixelType;
	  DoIt<InputPixelType, 2>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::INT:{
	  typedef int InputPixelType;
	  DoIt<InputPixelType, 2>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::ULONG:{
	  typedef unsigned long InputPixelType;
	  DoIt<InputPixelType, 2>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::LONG:{
	  typedef long InputPixelType;
	  DoIt<InputPixelType, 2>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::FLOAT:{
	  typedef float InputPixelType;
	  DoIt<InputPixelType, 2>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::DOUBLE:{
	  typedef double InputPixelType;
	  DoIt<InputPixelType, 2>(argc, argv);
	  break;
	}

	case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
	default:
	  std::cout << "unknown component type" << std::endl;
	  break;
	}//switch (componentType){
	break;
      
      case 3:
	switch (componentType){
	case itk::ImageIOBase::UCHAR:{
	  typedef unsigned char InputPixelType;
	  DoIt<InputPixelType, 3>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::CHAR:{
	  typedef char InputPixelType;
	  DoIt<InputPixelType, 3>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::USHORT:{
	  typedef unsigned short InputPixelType;
	  DoIt<InputPixelType, 3>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::SHORT:{
	  typedef short InputPixelType;
	  DoIt<InputPixelType, 3>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::UINT:{
	  typedef unsigned int InputPixelType;
	  DoIt<InputPixelType, 3>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::INT:{
	  typedef int InputPixelType;
	  DoIt<InputPixelType, 3>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::ULONG:{
	  typedef unsigned long InputPixelType;
	  DoIt<InputPixelType, 3>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::LONG:{
	  typedef long InputPixelType;
	  DoIt<InputPixelType, 3>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::FLOAT:{
	  typedef float InputPixelType;
	  DoIt<InputPixelType, 3>(argc, argv);
	  break;
	}
	case itk::ImageIOBase::DOUBLE:{
	  typedef double InputPixelType;
	  DoIt<InputPixelType, 3>(argc, argv);
	  break;
	}

	case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
	default:
	  std::cout << "unknown component type" << std::endl;
	  break;
	}//switch (componentType){
	break;
      
      default:
	std::cout << "Images of dimension " << dimensionType << " are not supported!" << std::endl;
	break;
      }//switch (dimensionType){

    }//try
    catch( itk::ExceptionObject &excep)
        {
        std::cerr << argv[0] << ": exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
        }
  
    return EXIT_SUCCESS;
    }

