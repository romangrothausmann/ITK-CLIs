////program to shift the labels of each slice such that the values do not overlap
//01: based on template.cxx
//02: try to use IterationEvent of SBS to access stat->GetMaximum()
//    could not figure out though how to store/add maximum of last slice for the next one


#include <complex>

#include "itkFilterWatcher.h" 
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkSliceBySliceImageFilter.h>
//#include <itkImageToImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkChangeLabelImageFilter.h>


unsigned long long m_lastMax= 0;


int dispatch_cT(itk::ImageIOBase::IOPixelType, itk::ImageIOBase::IOComponentType, size_t, int, char **);

template<typename InputComponentType>
int dispatch_pT(itk::ImageIOBase::IOPixelType pixelType, size_t, int, char **);

template<typename InputComponentType, typename InputPixelType>
int dispatch_D(size_t, int, char **);

template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int, char *argv[]);




template<typename InputImageType, typename OutputImageType, typename SBSInputImageType>
void FilterEventHandlerITK(itk::Object *caller, const itk::EventObject &event, void*){

    const itk::ProcessObject* filter = static_cast<const itk::ProcessObject*>(caller);

    if(itk::ProgressEvent().CheckEvent(&event))
	fprintf(stderr, "\r%s progress: %5.1f%%", filter->GetNameOfClass(), 100.0 * filter->GetProgress());//stderr is flushed directly
    else if(itk::IterationEvent().CheckEvent(&event)){
	std::cerr << " Iteration: " << (dynamic_cast<itk::SliceBySliceImageFilter<InputImageType, OutputImageType> *>(caller))->GetSliceIndex()
		  << " Max: " << +(dynamic_cast<itk::StatisticsImageFilter<SBSInputImageType> *>((dynamic_cast<itk::SliceBySliceImageFilter<InputImageType, OutputImageType> *>(caller))->GetInputFilter()))->GetMaximum() << std::endl;
	m_lastMax= (dynamic_cast<itk::StatisticsImageFilter<SBSInputImageType> *>((dynamic_cast<itk::SliceBySliceImageFilter<InputImageType, OutputImageType> *>(caller))->GetInputFilter()))->GetMaximum();
	}
    else if(itk::EndEvent().CheckEvent(&event))
	std::cerr << " Iteration: " << (dynamic_cast<itk::SliceBySliceImageFilter<InputImageType, OutputImageType> *>(caller))->GetSliceIndex()
		  << " Max: " << +(dynamic_cast<itk::StatisticsImageFilter<SBSInputImageType> *>((dynamic_cast<itk::SliceBySliceImageFilter<InputImageType, OutputImageType> *>(caller))->GetInputFilter()))->GetMaximum() << std::endl; 
    }



template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){


    typedef InputPixelType  OutputPixelType;
    
    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;


    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
 
    reader->SetFileName(argv[1]);
    // FilterWatcher watcherI(reader);
    // watcherI.QuietOn();
    // watcherI.ReportTimeOn();
    // try{ 
    //     reader->Update();
    //     }
    // catch(itk::ExceptionObject &ex){ 
    // 	std::cerr << ex << std::endl;
    // 	return EXIT_FAILURE;
    // 	}

    typename InputImageType::Pointer input= reader->GetOutput();



    typedef itk::SliceBySliceImageFilter<InputImageType, OutputImageType> SBSFilterType;
    typename SBSFilterType::Pointer sbs = SBSFilterType::New();
    sbs->SetInput(0, input);
    sbs->SetInput(1, input);

    typedef typename  SBSFilterType::InternalInputImageType SBSInputImageType;
    typedef typename  SBSFilterType::InternalOutputImageType SBSOutputImageType;


    // typedef itk::ImageToImageFilter<SBSInputImageType, SBSInputImageType> I2IType;
    // typename I2IType::Pointer I2I = I2IType::New();

    typedef itk::StatisticsImageFilter<SBSInputImageType> StatType;
    typename StatType::Pointer stat = StatType::New();

    typedef itk::AddImageFilter<SBSInputImageType, SBSInputImageType, SBSOutputImageType> AddType;
    typename AddType::Pointer adder = AddType::New();
    adder->SetInput1(stat->GetOutput());
    //adder->SetConstant2(stat->GetMaximum());
    adder->SetConstant2(m_lastMax);

    typedef itk::ChangeLabelImageFilter<SBSOutputImageType, SBSOutputImageType> ChangeLabType;
    typename ChangeLabType::Pointer ch= ChangeLabType::New();
    //ch->SetInput(adder->GetOutput());
    ch->SetInput(adder->GetOutput());
    ch->SetChange(m_lastMax, 0);

    //stat->SetInput(adder->GetOutput());

    sbs->SetInputFilter(stat);
    //sbs->SetInputFilter(adder);
    //sbs->SetOutputFilter(stat);
    sbs->SetOutputFilter(ch);

    itk::CStyleCommand::Pointer eventCallbackITK;
    eventCallbackITK = itk::CStyleCommand::New();
    eventCallbackITK->SetCallback(FilterEventHandlerITK<InputImageType, OutputImageType, SBSInputImageType>);

    sbs->AddObserver(itk::ProgressEvent(), eventCallbackITK);
    sbs->AddObserver(itk::IterationEvent(), eventCallbackITK);
    sbs->AddObserver(itk::EndEvent(), eventCallbackITK);
    try { 
	sbs->Update();
        }
    catch(itk::ExceptionObject &ex){ 
	std::cerr << ex << std::endl;
	return EXIT_FAILURE;
	}

    typename OutputImageType::Pointer output= sbs->GetOutput();

    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer);
    writer->SetFileName(argv[2]);
    writer->SetInput(output);
    //writer->UseCompressionOn();
    //writer->SetUseCompression(atoi(argv[3]));
    try{ 
        writer->Update();
        }
    catch(itk::ExceptionObject &ex){ 
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    return EXIT_SUCCESS;

    }


int dispatch_cT(itk::ImageIOBase::IOComponentType componentType, itk::ImageIOBase::IOPixelType pixelType, size_t dimensionType, int argc, char *argv[]){
  int res= 0;

  //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#a8dc783055a0af6f0a5a26cb080feb178
  //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00107
  //IOComponentType: UNKNOWNCOMPONENTTYPE, UCHAR, CHAR, USHORT, SHORT, UINT, INT, ULONG, LONG, FLOAT, DOUBLE

  switch (componentType){
  case itk::ImageIOBase::UCHAR:{
    typedef unsigned char InputComponentType;
    res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::CHAR:{
    typedef char InputComponentType;
    res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::USHORT:{
    typedef unsigned short InputComponentType;
    res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::SHORT:{
    typedef short InputComponentType;
    res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::UINT:{
    typedef unsigned int InputComponentType;
    res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::INT:{
    typedef int InputComponentType;
    res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::ULONG:{
    typedef unsigned long InputComponentType;
    res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::LONG:{
    typedef long InputComponentType;
    res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  // case itk::ImageIOBase::FLOAT:{
  //   typedef float InputComponentType;
  //   res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  // } break;
  // case itk::ImageIOBase::DOUBLE:{
  //   typedef double InputComponentType;
  //   res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  // } break;
  case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
  default:
    std::cout << "unknown component type" << std::endl;
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
  // case itk::ImageIOBase::RGB:{
  //   typedef itk::RGBPixel<InputComponentType> InputPixelType;
  //   res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
  // } break;
  // case itk::ImageIOBase::RGBA:{
  //   typedef itk::RGBAPixel<InputComponentType> InputPixelType;
  //   res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
  // } break;
  // case itk::ImageIOBase::COMPLEX:{
  //   typedef std::complex<InputComponentType> InputPixelType;
  //   res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
  // } break;
  // case itk::ImageIOBase::VECTOR:{
  //   typedef itk::VariableLengthVector<InputComponentType> InputPixelType;
  //   res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
  // } break;
  case itk::ImageIOBase::UNKNOWNPIXELTYPE:
  default:
    std::cerr << std::endl << "Error: Pixel type not handled!" << std::endl;
    break;
  }//switch 
  return res;
}


template<typename InputComponentType, typename InputPixelType>
int dispatch_D(size_t dimensionType, int argc, char *argv[]){
  int res= 0;
  switch (dimensionType){
  // case 1:
  //   res= DoIt<InputComponentType, InputPixelType, 1>(argc, argv);
  //   break;
  // case 2:
  //   res= DoIt<InputComponentType, InputPixelType, 2>(argc, argv);
  //   break;
  case 3:
    res= DoIt<InputComponentType, InputPixelType, 3>(argc, argv);
    break;
  default: 
    std::cerr << "Error: Images of dimension " << dimensionType << " are not handled!" << std::endl;
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
    typedef itk::Image<unsigned char, 3> ImageType;
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
    if ( argc != 3 ){
	std::cerr << "Missing Parameters: "
		  << argv[0]
		  << " Input_Image"
		  << " Output_Image"
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






