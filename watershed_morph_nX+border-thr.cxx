////program for iterative itkMorphologicalWatershedFromMarkersImageFilter
//01: based on watershed_morph_nX.cxx, added masking of borderImg

#include <string>
#include <sstream>

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
#include <itkGradientMagnitudeImageFilter.h>
#include <itkChangeLabelImageFilter.h>
#include <itkMaskImageFilter.h>



int dispatch_cT(itk::ImageIOBase::IOPixelType, itk::ImageIOBase::IOComponentType, size_t, int, char **);

template<typename InputComponentType>
int dispatch_pT(itk::ImageIOBase::IOPixelType pixelType, size_t, int, char **);

template<typename InputComponentType, typename InputPixelType>
int dispatch_D(size_t, int, char **);

template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int, char *argv[]);




template<typename ReaderImageType, typename WriterImageType>
void FilterEventHandlerITK(itk::Object *caller, const itk::EventObject &event, void*){

    const itk::ProcessObject* filter = static_cast<const itk::ProcessObject*>(caller);

    if(itk::ProgressEvent().CheckEvent(&event))
        fprintf(stderr, "\r%s progress: %5.1f%%", filter->GetNameOfClass(), 100.0 * filter->GetProgress());//stderr is flushed directly
    else if(itk::StartEvent().CheckEvent(&event)){
        if(strstr(filter->GetNameOfClass(), "ImageFileReader"))
            std::cerr << "Reading: " << (dynamic_cast<itk::ImageFileReader<ReaderImageType> *>(caller))->GetFileName();//cast only works if reader was instanciated for ReaderImageType!
        else if(strstr(filter->GetNameOfClass(), "ImageFileWriter"))
            std::cerr << "Writing: " << (dynamic_cast<itk::ImageFileWriter<WriterImageType> *>(caller))->GetFileName() << std::endl;//cast only works if writer was instanciated for WriterImageType!
        }
    else if(itk::EndEvent().CheckEvent(&event))
        std::cerr << std::endl;
    }


template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){

    typedef uint32_t  OutputPixelType;
    typedef uint8_t   MaskType;

    typedef itk::Image<InputPixelType, Dimension>   InputImageType;
    typedef itk::Image<double, Dimension>           GreyImageType;
    typedef itk::Image<OutputPixelType, Dimension>  LabelImageType;
    typedef itk::Image<MaskType,  Dimension>        MaskImageType;

    itk::CStyleCommand::Pointer eventCallbackITK;
    eventCallbackITK = itk::CStyleCommand::New();
    eventCallbackITK->SetCallback(FilterEventHandlerITK<InputImageType,LabelImageType>);


    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(argv[1]);
    reader->AddObserver(itk::AnyEvent(), eventCallbackITK);
    try{
        reader->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    typename GreyImageType::Pointer input;

    bool ws0_conn= true;//true reduces amount of watersheds
    bool ws_conn= false;

    double MinRelFacetSize= atof(argv[5]);
    uint8_t NumberOfExtraWS= atoi(argv[6]);
    char* interMedOutPrefix= NULL;

    if(argc == 9)
        interMedOutPrefix= argv[8];

    typename LabelImageType::Pointer markerImg;
    typename LabelImageType::Pointer borderImg;
    typename GreyImageType::Pointer gradientImg;
    typename LabelImageType::Pointer labelImg;
    typename LabelImageType::PixelType labelCnt;


    typedef itk::ImageFileWriter<LabelImageType>  LWriterType;
    typename LWriterType::Pointer lwriter = LWriterType::New();

    lwriter->SetUseCompression(atoi(argv[3]));
    lwriter->AddObserver(itk::AnyEvent(), eventCallbackITK);
    //lwriter->SetInput(labelImg);//doing that before each update to not get confused


    typedef itk::ShiftScaleImageFilter<InputImageType, GreyImageType> SSType;
    typename SSType::Pointer ss = SSType::New();
    if(atoi(argv[4]))
        ss->SetScale(-1); //invert by mul. with -1
    else
        ss->SetScale(1); //just convert to GreyImageType
    //ss->SetScale(atof(argv[4]));//abs(scaling) > 1 does not help gm to be more pronounced!
    ss->SetInput(reader->GetOutput());
    ss->Update();
    input= ss->GetOutput();

        {//scoped for better consistency
        typedef itk::HMinimaImageFilter<GreyImageType, GreyImageType> HMType; //seems for hmin in-type==out-type!!!
        typename HMType::Pointer hm= HMType::New();
        hm->SetHeight(MinRelFacetSize);
        hm->SetFullyConnected(ws0_conn);
        hm->SetInput(input);
        hm->AddObserver(itk::AnyEvent(), eventCallbackITK);
        hm->Update();

        typedef itk::RegionalMinimaImageFilter<GreyImageType, MaskImageType> RegMinType;
        typename RegMinType::Pointer rm = RegMinType::New();
        rm->SetFullyConnected(ws0_conn);
        rm->SetInput(hm->GetOutput());
        rm->AddObserver(itk::AnyEvent(), eventCallbackITK);
        rm->Update();

        // connected component labelling
        typedef itk::ConnectedComponentImageFilter<MaskImageType, LabelImageType> CCType;
        typename CCType::Pointer labeller = CCType::New();
        labeller->SetFullyConnected(ws0_conn);
        labeller->SetInput(rm->GetOutput());
        labeller->AddObserver(itk::AnyEvent(), eventCallbackITK);
        labeller->Update();
        labelImg= labeller->GetOutput();
        labelImg->DisconnectPipeline();
        labelCnt= labeller->GetObjectCount();
        }

    if(interMedOutPrefix){
        lwriter->SetFileName(std::string(interMedOutPrefix) + "_limg.mha");
        lwriter->SetInput(labelImg);
        lwriter->Update();
        }

    typedef itk::MorphologicalWatershedFromMarkersImageFilter<GreyImageType, LabelImageType> MWatershedType;
    typename MWatershedType::Pointer ws = MWatershedType::New();
    ws->SetMarkWatershedLine(NumberOfExtraWS); //use borders if higher order WS are wanted
    ws->SetFullyConnected(ws0_conn);
    ws->SetInput(input);
    ws->SetMarkerImage(labelImg);
    ws->AddObserver(itk::AnyEvent(), eventCallbackITK);
    ws->Update();

    if(NumberOfExtraWS > 0){
	if(interMedOutPrefix){
	    lwriter->SetFileName(std::string(interMedOutPrefix) + "_ws0.mha");
	    lwriter->SetInput(ws->GetOutput());
	    lwriter->Update();
	    }

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

	    typedef itk::BinaryThresholdImageFilter<InputImageType, MaskImageType> ThreshType2;
	    typename ThreshType2::Pointer th2 = ThreshType2::New();
	    th2->SetUpperThreshold(atof(argv[7]));
	    th2->SetInput(reader->GetOutput());
	    if(!atoi(argv[4])){ th2->SetOutsideValue(1); th2->SetInsideValue(0); }
	    else{ th2->SetOutsideValue(0); th2->SetInsideValue(1); }
	    th2->Update();

	    typedef itk::MaskImageFilter<LabelImageType, MaskImageType, LabelImageType> MaskType;
	    typename MaskType::Pointer mask = MaskType::New();
	    mask->SetInput(th->GetOutput());
	    mask->SetMaskImage(th2->GetOutput());
	    mask->ReleaseDataFlagOn();
	    mask->InPlaceOn();
	    mask->Update();

	    borderImg= mask->GetOutput();
	    borderImg->DisconnectPipeline();

	    if(interMedOutPrefix){
		lwriter->SetFileName(std::string(interMedOutPrefix) + "_bimg.mha");
		lwriter->SetInput(borderImg);
		lwriter->Update();
		}
	    }

	// to combine the markers again
	typedef itk::AddImageFilter<LabelImageType, LabelImageType, LabelImageType> AddType;
	typename AddType::Pointer adder = AddType::New();

	// to create gradient magnitude image
	typedef itk::GradientMagnitudeImageFilter<GreyImageType, GreyImageType> GMType;
	typename GMType::Pointer gm = GMType::New();

	gm->AddObserver(itk::AnyEvent(), eventCallbackITK);
	gradientImg= input;

	ws->SetMarkWatershedLine(false); //no use for a border in higher stages
	ws->SetFullyConnected(ws_conn);

	// to delete the background label
	typedef itk::ChangeLabelImageFilter<LabelImageType, LabelImageType> ChangeLabType;
	typename ChangeLabType::Pointer ch= ChangeLabType::New();
	ch->SetChange(labelCnt + 1, 0);


	for(char i= 0; i < NumberOfExtraWS; i++){

	    // Add the marker image to the watershed line image
	    adder->SetInput1(borderImg);
	    adder->SetInput2(labelImg);
	    adder->Update();
	    markerImg= adder->GetOutput();
	    markerImg->DisconnectPipeline();

	    // compute a gradient
	    gm->SetInput(gradientImg);
	    gm->Update();
	    gradientImg= gm->GetOutput();
	    gradientImg->DisconnectPipeline();

	    // Now apply higher order watershed
	    ws->SetInput(gradientImg);
	    ws->SetMarkerImage(markerImg);
	    ws->Update();

	    // delete the background label
	    ch->SetInput(ws->GetOutput());
	    ch->Update();
	    labelImg= ch->GetOutput();
	    labelImg->DisconnectPipeline();

	    if(interMedOutPrefix){
		std::stringstream sss;
		sss << interMedOutPrefix << "_ws" << i+1 << ".mha";
		lwriter->SetFileName(sss.str().c_str());
		//lwriter->SetFileName(std::string(interMedOutPrefix) + "_ws" + std::to_string(i+1) + ".mha");//c++11: std::to_string
		lwriter->SetInput(labelImg);
		lwriter->Update();
		}
	    }
        }
    else
        labelImg= ws->GetOutput();

    typedef itk::ImageFileWriter<LabelImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    writer->SetFileName(argv[2]);
    writer->SetInput(labelImg);
    writer->SetUseCompression(atoi(argv[3]));
    writer->AddObserver(itk::AnyEvent(), eventCallbackITK);
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
    if ( argc < 8 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Output_Image"
                  << " compress"
                  << " invert level WS-extra-runs"
                  << " threshold"
                  << " interMedOutPrefix"
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






