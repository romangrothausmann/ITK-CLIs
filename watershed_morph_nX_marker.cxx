////program for itkMorphologicalWatershedFromMarkersImageFilter
//01: based on template_2inputs.cxx and watershed_morph_nX_marker_UI8.cxx


#include <proc/readproc.h>//for look_up_our_self
#include <unistd.h>//for sysconf

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkCommand.h>

#include <itkShiftScaleImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkMorphologicalWatershedFromMarkersImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkChangeLabelImageFilter.h>



template<typename InputImageType>
void FilterEventHandlerITK(itk::Object *caller, const itk::EventObject &event, void*){

    const itk::ProcessObject* filter = static_cast<const itk::ProcessObject*>(caller);

    if(itk::ProgressEvent().CheckEvent(&event))
        fprintf(stderr, "\r%s progress: %5.1f%%", filter->GetNameOfClass(), 100.0 * filter->GetProgress());//stderr is flushed directly
    else if(itk::StartEvent().CheckEvent(&event)){
	if(strstr(filter->GetNameOfClass(), "ImageFileReader"))
	    std::cerr << "Reading: " << (dynamic_cast<itk::ImageFileReader<InputImageType> *>(caller))->GetFileName();
	}
    else if(itk::EndEvent().CheckEvent(&event))
        std::cerr << std::endl;
    }



template<typename InputComponent1, typename TypeInputComponentType2, typename InputPixelType1, typename InputPixelType2, size_t Dimension>
int DoIt(int argc, char *argv[]){

    typedef itk::Image<InputPixelType1, Dimension>  InputImageType1;
    typedef itk::Image<InputPixelType2, Dimension>  LabelImageType;
#ifdef USE_FLOAT
    typedef itk::Image<float, Dimension>            GreyImageType;
    std::cerr << "Using single precision (float)." << std::endl;
#else
    typedef itk::Image<double, Dimension>           GreyImageType;
    std::cerr << "Using double precision (double)." << std::endl;
#endif

    itk::CStyleCommand::Pointer eventCallbackITK;
    eventCallbackITK = itk::CStyleCommand::New();
    eventCallbackITK->SetCallback(FilterEventHandlerITK<InputImageType1>);//only works for both readers if InputImageType1 == InputImageType2

    ////for mem monitoring: http://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-run-time-in-c
    struct proc_t usage;//description in: /usr/include/proc/readproc.h
    double page_size_mb = sysconf(_SC_PAGE_SIZE) / 1024. / 1024.; // in case x86-64 is configured to use 2MB pages


    typename GreyImageType::Pointer input;
        {
        typedef itk::ImageFileReader<InputImageType1> ReaderType1;
        typename ReaderType1::Pointer reader1 = ReaderType1::New();

        reader1->SetFileName(argv[1]);
        reader1->AddObserver(itk::AnyEvent(), eventCallbackITK);
        try{
            reader1->Update();
            }
        catch(itk::ExceptionObject &ex){
            std::cerr << ex << std::endl;
            return EXIT_FAILURE;
            }

        typedef itk::ShiftScaleImageFilter<InputImageType1, GreyImageType> SSType;
        typename SSType::Pointer ss = SSType::New();
        if(atoi(argv[6]))
            ss->SetScale(-1); //invert by mul. with -1
        else
            ss->SetScale(1); //just convert to GreyImageType
        ss->SetInput(reader1->GetOutput());
        ss->AddObserver(itk::ProgressEvent(), eventCallbackITK);
        ss->AddObserver(itk::EndEvent(), eventCallbackITK);
        try{
            ss->Update();
            }
        catch(itk::ExceptionObject &ex){
            std::cerr << ex << std::endl;
            return EXIT_FAILURE;
            }
        input= ss->GetOutput();
        input->DisconnectPipeline();//will need its own Delete later on!
        }
    look_up_our_self(&usage); fprintf(stderr, "vsize: %.3f mb; rss: %.3f mb\n", usage.vsize/1024./1024., usage.rss * page_size_mb);

    typename LabelImageType::Pointer labelImg;
    typename LabelImageType::PixelType labelCnt;
        {
        typedef itk::ImageFileReader<LabelImageType> ReaderType2;
        typename ReaderType2::Pointer reader2 = ReaderType2::New();

        reader2->SetFileName(argv[2]);
        //reader2->AddObserver(itk::AnyEvent(), eventCallbackITK);
        try{
            reader2->Update();
            }
        catch(itk::ExceptionObject &ex){
            std::cerr << ex << std::endl;
            return EXIT_FAILURE;
            }

        typedef itk::StatisticsImageFilter<LabelImageType> FilterType;
        typename FilterType::Pointer stat= FilterType::New();
        stat->SetInput(reader2->GetOutput());
        //stat->InPlaceOn();//not available
        stat->ReleaseDataFlagOn();

        stat->AddObserver(itk::ProgressEvent(), eventCallbackITK);
        stat->AddObserver(itk::EndEvent(), eventCallbackITK);
        try{
            stat->Update();
            }
        catch(itk::ExceptionObject &ex){
            std::cerr << ex << std::endl;
            return EXIT_FAILURE;
            }

        std::cerr << "Min: " << +stat->GetMinimum() << " Max: " << +stat->GetMaximum() << " Mean: " << +stat->GetMean() << " Std: " << +stat->GetSigma() << " Variance: " << +stat->GetVariance() << " Sum: " << +stat->GetSum() << std::endl;

        labelCnt= stat->GetMaximum();
        labelImg= stat->GetOutput();
        labelImg->DisconnectPipeline();//will need its own Delete later on!
        }
    look_up_our_self(&usage); fprintf(stderr, "vsize: %.3f mb; rss: %.3f mb\n", usage.vsize/1024./1024., usage.rss * page_size_mb);

    bool ws0_conn= true;//true reduces amount of watersheds
    bool ws_conn= ws0_conn;
    if(argc == 9)
	ws_conn= atoi(argv[8]);

    uint8_t NumberOfExtraWS= atoi(argv[5]);

    typename LabelImageType::Pointer markerImg;
    typename LabelImageType::Pointer borderImg;
    typename GreyImageType::Pointer gradientImg;


    typedef itk::MorphologicalWatershedFromMarkersImageFilter<GreyImageType, LabelImageType> MWatershedType;
    typename MWatershedType::Pointer ws = MWatershedType::New();
    ws->SetMarkWatershedLine(NumberOfExtraWS); //use borders if higher order WS are wanted
    ws->SetFullyConnected(ws0_conn);
    ws->SetInput(input);
    ws->SetMarkerImage(labelImg);
    ws->AddObserver(itk::ProgressEvent(), eventCallbackITK);
    ws->AddObserver(itk::EndEvent(), eventCallbackITK);
    ws->Update();

    look_up_our_self(&usage); fprintf(stderr, "vsize: %.3f mb; rss: %.3f mb\n", usage.vsize/1024./1024., usage.rss * page_size_mb);
    if(NumberOfExtraWS > 0){

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

        look_up_our_self(&usage); fprintf(stderr, "vsize: %.3f mb; rss: %.3f mb\n", usage.vsize/1024./1024., usage.rss * page_size_mb);
        // to combine the markers again
        typedef itk::AddImageFilter<LabelImageType, LabelImageType, LabelImageType> AddType;
        typename AddType::Pointer adder = AddType::New();
        adder->InPlaceOn();
        adder->ReleaseDataFlagOn();

        // to create gradient magnitude image
        typedef itk::GradientMagnitudeImageFilter<GreyImageType, GreyImageType> GMType;
        typename GMType::Pointer gm = GMType::New();

        gm->AddObserver(itk::ProgressEvent(), eventCallbackITK);
        gm->AddObserver(itk::EndEvent(), eventCallbackITK);
        //gm->InPlaceOn();//not available
        //gm->ReleaseDataFlagOn();//gm output is used for ws and next gm!
        gradientImg= input;
        input->ReleaseDataFlagOn();//free input as soon as gradientImg has been used as input by gm

        ws->SetMarkWatershedLine(atoi(argv[7]));//border in higher stages can make a difference but will not separate lables that initially touch already!
        ws->SetFullyConnected(ws_conn);
        //ws->InPlaceOn();//not available
        //ws->ReleaseDataFlagOn();//no problem but not needed as ch->InPlaceOn()

        // to delete the background label
        typedef itk::ChangeLabelImageFilter<LabelImageType, LabelImageType> ChangeLabType;
        typename ChangeLabType::Pointer ch= ChangeLabType::New();
        ch->SetChange(labelCnt + 1, 0);
        ch->InPlaceOn();
        ch->ReleaseDataFlagOn();//will be handled by adder if adder->InPlaceOn() AND used for adder->SetInput1

        std::cerr << "Starting extra runs..." << std::endl;

        for(char i= 0; i < NumberOfExtraWS; i++){
	    std::cerr << "Run: " << i+1 << std::endl;
            //// DisconnectPipeline() on outputs does not make sense here becauses:
            //// - the filter it belongs to is not scoped to the loop
            //// - it likely avoids mem freeing expected from ReleaseDataFlagOn()

            // Add the marker image to the watershed line image
            //// labelImg will not be needed again afterwards
            //// so setting Input1 to labelImg with InPlaceOn() will overwrite labelImg
            //// with InPlaceOn() use Input2 for borderImg to avoid loosing orig borderImg
            adder->SetInput1(labelImg);//if InPlaceOn(), input1 will be changed! and used as output!
            adder->SetInput2(borderImg);
            adder->Update();//frees mem of labelImg if ch->ReleaseDataFlagOn(); even if adder->InPlaceOn();?
            markerImg= adder->GetOutput();
            //markerImg->DisconnectPipeline();
            //labelImg->Delete();//free mem of labelImg originating from initial markers; do not use if adder->InPlaceOn()

            look_up_our_self(&usage); fprintf(stderr, "vsize: %.3f mb; rss: %.3f mb\n", usage.vsize/1024./1024., usage.rss * page_size_mb);
	    gradientImg->ReleaseDataFlagOn();
            // compute a gradient
            gm->SetInput(gradientImg);
            gm->Update();
            //gradientImg->Delete();//free mem of orig input / last gradientImg <- bad! causes double free! use ReleaseDataFlag on reader instead and rely on smart ponter logic for last gm-output? http://public.kitware.com/pipermail/insight-users/2009-October/033004.html
            gradientImg= gm->GetOutput();
            gradientImg->DisconnectPipeline();//segfaults without! Why?

            look_up_our_self(&usage); fprintf(stderr, "vsize: %.3f mb; rss: %.3f mb\n", usage.vsize/1024./1024., usage.rss * page_size_mb);
	    markerImg->ReleaseDataFlagOn();
            // Now apply higher order watershed
            ws->SetInput(gradientImg);
            ws->SetMarkerImage(markerImg);
            ws->Update();//frees mem of markerImg if adder->ReleaseDataFlagOn();

            look_up_our_self(&usage); fprintf(stderr, "vsize: %.3f mb; rss: %.3f mb\n", usage.vsize/1024./1024., usage.rss * page_size_mb);
            // delete the background label
            ch->SetInput(ws->GetOutput());//with ch->InPlaceOn() ws output will be overwritten!
            ch->Update();//frees mem of ws output if ws->ReleaseDataFlagOn();
            labelImg= ch->GetOutput();
            //labelImg->DisconnectPipeline();
            }
        }
    else
        labelImg= ws->GetOutput();

    typedef itk::ImageFileWriter<LabelImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    writer->SetFileName(argv[3]);
    writer->SetInput(labelImg);
    writer->SetUseCompression(atoi(argv[4]));
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


template<typename InputComponentType1, typename InputComponentType2, typename InputPixelType1, typename InputPixelType2>
int dispatch_D(size_t dimensionType, int argc, char *argv[]){
    int res= 0;
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
    int res= 0;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (pixelType2){
    case itk::ImageIOBase::SCALAR:{
        typedef InputComponentType2 InputPixelType2;
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
    int res= 0;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (pixelType1){
    case itk::ImageIOBase::SCALAR:{
        typedef InputComponentType1 InputPixelType1;
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
    int res= 0;

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
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
    default:
        std::cerr << "unknown component type" << std::endl;
        break;
        }//switch
    return res;
    }

int dispatch_cT1(itk::ImageIOBase::IOComponentType componentType1, itk::ImageIOBase::IOComponentType componentType2, itk::ImageIOBase::IOPixelType pixelType1, itk::ImageIOBase::IOPixelType pixelType2, size_t dimensionType, int argc, char *argv[]){
    int res= 0;

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
    if ( argc < 8 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Marker_Image"
                  << " Output_Image"
                  << " compress"
                  << " NumberOfExtraWS invert"
                  << " borderExtraWS"
                  << " [connExtraWS]"
                  << std::endl;

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






