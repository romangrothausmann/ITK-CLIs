////program for iterative itkMorphologicalWatershedFromMarkersImageFilter
//01: based on watershed_morph_nX, but grows borders on HMin eliminated regions (does not use gradient image)


#include <complex>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkCommand.h>

#include <itkShiftScaleImageFilter.h>
#include <itkHMinimaImageFilter.h>
#include <itkRegionalMinimaImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkMorphologicalWatershedFromMarkersImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkChangeLabelImageFilter.h>



int dispatch_cT(itk::ImageIOBase::IOPixelType, itk::ImageIOBase::IOComponentType, size_t, int, char **);

template<typename InputComponentType>
int dispatch_pT(itk::ImageIOBase::IOPixelType pixelType, size_t, int, char **);

template<typename InputComponentType, typename InputPixelType>
int dispatch_D(size_t, int, char **);

template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int, char *argv[]);




void FilterEventHandlerITK(itk::Object *caller, const itk::EventObject &event, void*){

    const itk::ProcessObject* filter = static_cast<const itk::ProcessObject*>(caller);

    if(itk::ProgressEvent().CheckEvent(&event))
	fprintf(stderr, "\r%s progress: %5.1f%%", filter->GetNameOfClass(), 100.0 * filter->GetProgress());//stderr is flushed directly
    else if(itk::EndEvent().CheckEvent(&event))
	std::cerr << std::endl << std::flush;   
    }


template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){

    typedef uint32_t  OutputPixelType;
    typedef uint8_t   MaskType;

    typedef itk::Image<InputPixelType, Dimension>   GreyImageType;
    typedef itk::Image<OutputPixelType, Dimension>  LabelImageType;
    typedef itk::Image<MaskType,  Dimension>        MaskImageType;

    itk::CStyleCommand::Pointer eventCallbackITK;
    eventCallbackITK = itk::CStyleCommand::New();
    eventCallbackITK->SetCallback(FilterEventHandlerITK);


    typedef itk::ImageFileReader<GreyImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
 
    reader->SetFileName(argv[1]);
    reader->AddObserver(itk::ProgressEvent(), eventCallbackITK);
    reader->AddObserver(itk::EndEvent(), eventCallbackITK);
    try{ 
        reader->Update();
        }
    catch(itk::ExceptionObject &ex){ 
	std::cerr << ex << std::endl;
	return EXIT_FAILURE;
	}

    typename GreyImageType::Pointer input= reader->GetOutput();

    bool ws0_conn= true;//true reduces amount of watersheds
    bool ws_conn= false;

    double MinRelFacetSize= atof(argv[3]);

    typename LabelImageType::Pointer markerImg;
    typename LabelImageType::Pointer borderImg;
    typename LabelImageType::Pointer labelImg;
    typename LabelImageType::PixelType labelCnt;

    // typedef itk::ShiftScaleImageFilter<GreyImageType, GreyImageType> SSType;
    // SSType::Pointer ss = SSType::New();
    // ss->SetScale(-1); //invert by mul. with -1
    // ss->SetInput(input);
    // ss->Update();
    // input= ss->GetOutput();

    {//scoped for better consistency
    typedef itk::HMinimaImageFilter<GreyImageType, GreyImageType> HMType; //seems for hmin in-type==out-type!!!
    typename HMType::Pointer hm= HMType::New();
    hm->SetHeight(MinRelFacetSize);
    hm->SetFullyConnected(ws0_conn);
    hm->SetInput(input);
    hm->AddObserver(itk::ProgressEvent(), eventCallbackITK);
    hm->AddObserver(itk::EndEvent(), eventCallbackITK);
    hm->Update();

    typedef itk::RegionalMinimaImageFilter<GreyImageType, MaskImageType> RegMinType;
    typename RegMinType::Pointer rm = RegMinType::New();
    rm->SetFullyConnected(ws0_conn);
    rm->SetInput(hm->GetOutput());
    rm->AddObserver(itk::ProgressEvent(), eventCallbackITK);
    rm->AddObserver(itk::EndEvent(), eventCallbackITK);
    rm->Update();

    // connected component labelling
    typedef itk::ConnectedComponentImageFilter<MaskImageType, LabelImageType> CCType;
    typename CCType::Pointer labeller = CCType::New();
    labeller->SetFullyConnected(ws0_conn);
    labeller->SetInput(rm->GetOutput());
    labeller->AddObserver(itk::ProgressEvent(), eventCallbackITK);
    labeller->AddObserver(itk::EndEvent(), eventCallbackITK);
    labeller->Update();
    labelImg= labeller->GetOutput();
    labelImg->DisconnectPipeline();
    labelCnt= labeller->GetObjectCount();
    }

    typedef itk::MorphologicalWatershedFromMarkersImageFilter<GreyImageType, LabelImageType> MWatershedType;
    typename MWatershedType::Pointer ws = MWatershedType::New();
    ws->SetMarkWatershedLine(true); //use borders if higher order WS are wanted
    ws->SetFullyConnected(ws0_conn);
    ws->SetInput(input);
    ws->SetMarkerImage(labelImg);
    ws->AddObserver(itk::ProgressEvent(), eventCallbackITK);
    ws->AddObserver(itk::EndEvent(), eventCallbackITK);
    ws->Update(); 

    {//scoped for better consistency
    // extract the watershed lines and combine with the orginal markers
    typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType> ThreshType;
    typename ThreshType::Pointer th = ThreshType::New();
    th->SetUpperThreshold(0);
    th->SetOutsideValue(0);
    // set the inside value to the number of markers + 1
    th->SetInsideValue(labelCnt + 1);
    th->SetInput(ws->GetOutput());
    th->Update();
    borderImg= th->GetOutput();
    borderImg->DisconnectPipeline();
    }

    // to combine the markers again
    typedef itk::AddImageFilter<LabelImageType, LabelImageType, LabelImageType> AddType;
    typename AddType::Pointer adder = AddType::New();

    ws->SetMarkWatershedLine(false); //no use for a border in higher stages
    ws->SetFullyConnected(ws_conn); 

    // to delete the background label
    typedef itk::ChangeLabelImageFilter<LabelImageType, LabelImageType> ChangeLabType;
    typename ChangeLabType::Pointer ch= ChangeLabType::New();
    ch->SetChange(labelCnt + 1, 0);

    // Add the marker image to the watershed line image
    adder->SetInput1(borderImg);
    adder->SetInput2(labelImg);
    adder->Update();
    markerImg= adder->GetOutput();
    markerImg->DisconnectPipeline();
    
    // Now apply higher order watershed
    ws->SetMarkerImage(markerImg);
    ws->Update();
    
    // delete the background label
    ch->SetInput(ws->GetOutput());
    ch->Update();
    labelImg= ch->GetOutput();
    labelImg->DisconnectPipeline();

    typedef itk::ImageFileWriter<LabelImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    writer->SetFileName(argv[2]);
    writer->SetInput(labelImg);
    writer->UseCompressionOn();
    //writer->SetUseCompression(atoi(argv[]));
    writer->AddObserver(itk::ProgressEvent(), eventCallbackITK);
    writer->AddObserver(itk::EndEvent(), eventCallbackITK);
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
  case itk::ImageIOBase::FLOAT:{
    typedef float InputComponentType;
    res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
  case itk::ImageIOBase::DOUBLE:{
    typedef double InputComponentType;
    res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
  } break;
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
		  << " MinRelFacetSize"
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






