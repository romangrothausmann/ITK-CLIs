////program to shift the labels of each slice such that the values do not overlap
//01: based on Modules/Filtering/ImageCompose/test/itkJoinSeriesImageFilterStreamingTest.cxx


#include <complex>

#include "itkFilterWatcher.h" 
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkExtractImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkChangeLabelImageFilter.h>
#include <itkJoinSeriesImageFilter.h>
#include <itkPipelineMonitorImageFilter.h>


int dispatch_cT(itk::ImageIOBase::IOPixelType, itk::ImageIOBase::IOComponentType, size_t, int, char **);

template<typename InputComponentType>
int dispatch_pT(itk::ImageIOBase::IOPixelType pixelType, size_t, int, char **);

template<typename InputComponentType, typename InputPixelType>
int dispatch_D(size_t, int, char **);

template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int, char *argv[]);




// template<typename InputImageType, typename OutputImageType>
// void FilterEventHandlerITK(itk::Object *caller, const itk::EventObject &event, void*){

//     const itk::ProcessObject* filter = static_cast<const itk::ProcessObject*>(caller);

//     if(itk::ProgressEvent().CheckEvent(&event))
// 	fprintf(stderr, "\r%s progress: %5.1f%%", filter->GetNameOfClass(), 100.0 * filter->GetProgress());//stderr is flushed directly
//     else if(itk::IterationEvent().CheckEvent(&event))
//      std::cerr << " Iteration: " << (dynamic_cast<itk::SliceBySliceImageFilter<InputImageType, OutputImageType> *>(caller))->GetSliceIndex() << std::endl;   
//     else if(strstr(filter->GetNameOfClass(), "ImageFileReader"))
// 	std::cerr << "Reading: " << (dynamic_cast<itk::ImageFileReader<InputImageType> *>(caller))->GetFileName() << std::endl;   
//     else if(itk::EndEvent().CheckEvent(&event))
// 	std::cerr << std::endl;   
//     }



template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){

    typedef uint32_t  OutputPixelType;
    
    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;

    typedef itk::Image<OutputPixelType, Dimension-1> SliceImageType;

    // itk::CStyleCommand::Pointer eventCallbackITK;
    // eventCallbackITK = itk::CStyleCommand::New();
    // eventCallbackITK->SetCallback(FilterEventHandlerITK<InputImageType, OutputImageType>);


    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
 
    reader->SetFileName(argv[1]);

    reader->UpdateOutputInformation();
    typename InputImageType::RegionType region= reader->GetOutput()->GetLargestPossibleRegion();

    std::cerr << "input region index: " << region.GetIndex()
	      << "  size: " <<  region.GetSize()
	      << std::endl;

    
    typedef itk::ExtractImageFilter<InputImageType,SliceImageType>     SliceExtractorFilterType;
    typedef itk::StatisticsImageFilter<SliceImageType> StatType;
    typedef itk::AddImageFilter<SliceImageType, SliceImageType> AddType;
    typedef itk::ChangeLabelImageFilter<SliceImageType, SliceImageType> ChangeLabType;
    typedef itk::JoinSeriesImageFilter<SliceImageType, OutputImageType> JoinSeriesFilterType;
    typedef itk::ImageFileWriter<OutputImageType>                       ImageFileWriterType;

    const unsigned int numberOfSlices = itk::Math::CastWithRangeCheck<unsigned int>(reader->GetOutput()->GetLargestPossibleRegion().GetSize(Dimension-1));


    typename itk::PipelineMonitorImageFilter<InputImageType>::Pointer monitor1 = itk::PipelineMonitorImageFilter<InputImageType>::New();
    monitor1->SetInput( reader->GetOutput() );

    std::vector<itk::ProcessObject::Pointer> savedPointers;

    typename JoinSeriesFilterType::Pointer joinSeries = JoinSeriesFilterType::New();
    joinSeries->SetOrigin( reader->GetOutput()->GetOrigin()[Dimension-1] );
    joinSeries->SetSpacing( reader->GetOutput()->GetSpacing()[Dimension-1] );

    OutputPixelType m_lastMax= 0;

    for (typename InputImageType::SizeValueType z = 0; z < numberOfSlices; ++z )
	{

	typename SliceExtractorFilterType::Pointer extractor = SliceExtractorFilterType::New();
	extractor->SetDirectionCollapseToSubmatrix();

	typename SliceExtractorFilterType::InputImageRegionType slice( reader->GetOutput()->GetLargestPossibleRegion() );
	slice.SetSize( Dimension-1, 0 );
	slice.SetIndex( Dimension-1, z );

	std::cerr << "slice region index: " << slice.GetIndex()
		  << "  size: " << slice.GetSize()
		  << std::endl;


	extractor->SetExtractionRegion( slice );
	extractor->SetInput( monitor1->GetOutput() );
	extractor->InPlaceOn();
	extractor->ReleaseDataFlagOn();

	typename StatType::Pointer stat = StatType::New();
	stat->SetInput(extractor->GetOutput());
	//FilterWatcher watcher2(stat);

	typename AddType::Pointer adder = AddType::New();
	adder->SetInput1(stat->GetOutput());
	adder->SetConstant2(m_lastMax);
	//FilterWatcher watcher3(adder);

	typename ChangeLabType::Pointer ch= ChangeLabType::New();
	ch->SetInput(adder->GetOutput());
	ch->SetChange(m_lastMax, 0); //shift BG back to 0
	//FilterWatcher watcher4(ch);

	ch->Update();
	//stat->Update();

	savedPointers.push_back( ch.GetPointer() );

	joinSeries->PushBackInput( ch->GetOutput() );

	m_lastMax+= stat->GetMaximum();
	std::cerr << "m_lastMax: " << +m_lastMax << std::endl;

	}


    typename itk::PipelineMonitorImageFilter<OutputImageType>::Pointer monitor2 = itk::PipelineMonitorImageFilter<OutputImageType>::New();
    monitor2->SetInput( joinSeries->GetOutput() );

    typename ImageFileWriterType::Pointer writer = ImageFileWriterType::New();
    writer->SetInput( monitor2->GetOutput() );
    writer->SetFileName(argv[2]);
    writer->SetUseCompression(atoi(argv[3]));
    writer->SetNumberOfStreamDivisions( numberOfSlices );


    try
	{
	writer->Update();
	}
    catch (...)
	{
	std::cerr << "Exception while trying to stream write file." << std::endl;
	throw;
	}

    std::cout << "Number of Updates: " << monitor1->GetNumberOfUpdates() << std::endl;
    std::cout << "Verifying ImageFileReader to ExtractImageFilter pipeline interaction" << std::endl;

    // We can not use one of the standard verify all methods due to
    // multiple filters connected to the output of the reader
    if ( !(monitor1->VerifyInputFilterExecutedStreaming( numberOfSlices ) &&
	    monitor1->VerifyInputFilterMatchedUpdateOutputInformation()) )
	{
	//std::cerr << monitor1;
	//return EXIT_FAILURE;
	}

    std::cout << "Verifying JoinSeriesImageFilter to ImageFileWriter pipeline interaction" << std::endl;
    if ( !monitor2->VerifyAllInputCanStream( numberOfSlices ) )
	{
	std::cerr << monitor2;
	//return EXIT_FAILURE;
	}

    return EXIT_SUCCESS;

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
    if ( argc != 4 ){
	std::cerr << "Missing Parameters: "
		  << argv[0]
		  << " Input_Image"
		  << " Output_Image"
		  << " compress"
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






