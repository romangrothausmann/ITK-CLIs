#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkShiftScaleImageFilter.h>
#include <itkMorphologicalWatershedFromMarkersImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkFlatStructuringElement.h>
#include <itkGrayscaleErodeImageFilter.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkBlackTopHatImageFilter.h>
#include <itkBinaryImageToShapeLabelMapFilter.h>
#include <itkShapeOpeningLabelMapFilter.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkParabolicErodeImageFilter.h>
#include <itkLabelMapToLabelImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkMaximumImageFilter.h>
#include <itkLabelStatisticsImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkChangeLabelImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include "tclap/CmdLine.h"

bool WriteDebug = false;
bool Compression = false;
std::string DebugDir("./");

template <class TImage>
void writeIm(typename TImage::Pointer Im, std::string filename){
    typedef typename itk::ImageFileWriter<TImage> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    FilterWatcher watcherO(writer);
    writer->SetInput(Im);
    writer->SetUseCompression(Compression);
    writer->SetFileName(filename.c_str());
    writer->Update();
    }


template <class ImType>
void writeImDbg(typename ImType::Pointer Im, std::string filename){
    if (WriteDebug){
	writeIm<ImType>(Im, DebugDir + filename);
	}
    }


typedef class CmdLineType{
public:
    std::string InputImFile, OutputImFile;
  int DarkVol, markererode, markervolume, MorphGradRad, darkmarkerdilate;
    float SmoothGradSigma;
    bool Compression;
    } CmdLineType;

void ParseCmdLine(int argc, char* argv[], CmdLineType &CmdLineObj){
    using namespace TCLAP;
    try {
	// Define the command line object.
	CmdLine cmd("Lung segmentation prototype ", ' ', "0.9");


	ValueArg<std::string> inArg("i","input","input image",true,"result","string");
	cmd.add( inArg );

	ValueArg<std::string> outArg("o","output","output image", true,"","string");
	cmd.add( outArg );

	// ValueArg<int> darkArg("","darkvol","Volume used for attribute closing (voxels)",true, 3000, "integer");
	// cmd.add( darkArg);

	ValueArg<int> markererodeArg("","markererode","Radius of erosion of thresholded image used to produce markers (voxels)",true, 30, "integer");
	cmd.add( markererodeArg);

	ValueArg<int> darkmarkerdilateArg("","darkmarkerdilate","Radius of dilation of image used to produce bloodcell markers (voxels)",true, 20, "integer");
	cmd.add( darkmarkerdilateArg );

	ValueArg<int> markervolArg("","markervol","Minimum volume of a bright marker",false, 100000, "integer");
	cmd.add( markervolArg);

	ValueArg<int> morphgradArg("","morphgradrad","Radius of SE for morphologicad gradient",false, 1, "integer");
	cmd.add( morphgradArg );

	ValueArg<float> gradsmoothArg("","gradsigma","Smoothing kernel (voxels)",false, 1, "float");
	cmd.add( gradsmoothArg);

	SwitchArg compArg("c", "compress", "write compressed images", false);
	cmd.add( compArg );

	SwitchArg dbgArg("d", "debug", "write debug images", false);
	cmd.add( dbgArg );

	ValueArg<std::string> dfArg("","debugfolder","location for debug images", false,"/tmp/","string");
	cmd.add( dfArg );


	// Parse the args.
	cmd.parse( argc, argv );

	// Get the value parsed by each arg.
	CmdLineObj.InputImFile = inArg.getValue();
	CmdLineObj.OutputImFile = outArg.getValue();
	//CmdLineObj.DarkVol = darkArg.getValue();
	CmdLineObj.Compression = compArg.getValue();
	CmdLineObj.markererode = markererodeArg.getValue();
	CmdLineObj.markervolume = markervolArg.getValue();
	CmdLineObj.MorphGradRad = morphgradArg.getValue();
	CmdLineObj.SmoothGradSigma = gradsmoothArg.getValue();
	CmdLineObj.darkmarkerdilate = darkmarkerdilateArg.getValue();
	Compression = CmdLineObj.Compression;
	WriteDebug = dbgArg.getValue();
	DebugDir = dfArg.getValue();
	}
    catch (ArgException &e)  // catch any exceptions
	{
	std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	}
    }

template <class InImType, class OutputImType> 
typename OutputImType::Pointer findBrightMarkers(
    typename InImType::Pointer input, 
    int erad,
    int volumeopening,
    int &MaxLabel){
    // basic marker finding. Threshold, erode binary and keep larger
    // blobs
    typedef typename itk::Image<unsigned char, InImType::ImageDimension > BinaryType;
    typedef typename itk::OtsuThresholdImageFilter<InImType, BinaryType> OtsuType;
    typename OtsuType::Pointer Otsu = OtsuType::New();
    Otsu->SetInput(input);
    Otsu->SetInsideValue(0);
    Otsu->SetOutsideValue(1);

    typedef typename itk::ParabolicErodeImageFilter<BinaryType, BinaryType> EType;
    typename EType::Pointer eroder = EType::New();
    eroder->SetInput(Otsu->GetOutput());
    eroder->SetUseImageSpacing(false);
    eroder->SetScale(erad);

    typedef typename itk::BinaryImageToShapeLabelMapFilter<BinaryType> LabellerType;
    typename LabellerType::Pointer labeller = LabellerType::New();
    labeller->SetInput(eroder->GetOutput());
    labeller->SetInputForegroundValue(1);
    typedef typename itk::ShapeOpeningLabelMapFilter<typename LabellerType::OutputImageType> ShapeFilterType;
    typename ShapeFilterType::Pointer shapefilter = ShapeFilterType::New();
    shapefilter->SetLambda(volumeopening);
    shapefilter->SetAttribute("NumberOfPixels");
    shapefilter->SetInput(labeller->GetOutput());
  
    typedef typename itk::LabelMapToLabelImageFilter<typename LabellerType::OutputImageType, OutputImType> ConvType;

    typename ConvType::Pointer tolabIm = ConvType::New();
    tolabIm->SetInput(shapefilter->GetOutput());
  
    typename OutputImType::Pointer result = tolabIm->GetOutput();
    result->Update();
    result->DisconnectPipeline();
    MaxLabel = labeller->GetOutput()->GetNumberOfLabelObjects();
    return(result);
    }

template <class InImType, class OutputImType> 
typename OutputImType::Pointer findDarkMarkersWS(
    typename InImType::Pointer input,
    typename OutputImType::Pointer markers,
    int BorderLabel){
    // Tesselate the bright regions by doing a watershed on the inverted
    // input.
    // Invert by subtracting input from the max, so no funny business
    // with unsigned types.

    typedef typename itk::StatisticsImageFilter<InImType> StatsType;
    typename StatsType::Pointer stats = StatsType::New();
    stats->SetInput(input);
    stats->Update();

    typedef typename itk::ShiftScaleImageFilter<InImType, InImType> ShiftType;
    typename ShiftType::Pointer shifter = ShiftType::New();
    shifter->SetInput(input);
    shifter->SetShift(- stats->GetMaximum() );
    shifter->SetScale(-1);

    typedef typename itk::MorphologicalWatershedFromMarkersImageFilter<InImType, OutputImType> WSType;
    typename WSType::Pointer ws = WSType::New();
    FilterWatcher watcherws1(ws);
    ws->SetInput(shifter->GetOutput());
    ws->SetMarkerImage(markers);
    ws->SetMarkWatershedLine(true);

    typedef typename itk::BinaryThresholdImageFilter<OutputImType, OutputImType> SelectType;
    typename SelectType::Pointer select = SelectType::New();
    select->SetInput(ws->GetOutput());
    select->SetLowerThreshold(0);
    select->SetUpperThreshold(0);
    select->SetInsideValue(BorderLabel);
    select->SetOutsideValue(0);
    // Sometimes a foreground marker will be incorrectly broken in
    // two, leading to an undesired boundary splitting them after the
    // first stage watershed. This will make it difficult to join
    // regions in the interactive phase. Can try removing pixels from
    // the dark markers, if they are too bright. The importance of this
    // step will depend on a combination of the erosion size and topology.
    // Take mean and SD of markers, discard borders that are bright by
    // this metric.

    // compute label stats of the marker image - should be a binary, but
    // I can't be bothered thresholding again. We'll just use the stats
    // of the biggest marker. TODO combine mean, count and SD - look up
    // the formula
    typedef typename itk::LabelStatisticsImageFilter<InImType, OutputImType> LabStatsType;

    typename LabStatsType::Pointer markstats = LabStatsType::New();
    markstats->SetInput(input);
    markstats->SetLabelInput(markers);
    markstats->Update();
  
    typedef typename LabStatsType::ValidLabelValuesContainerType ValidLabelValuesType;
    int msize = 0;
    typename LabStatsType::RealType MMean(0), MSigma(0);

    for (typename ValidLabelValuesType::const_iterator vIt=markstats->GetValidLabelValues().begin();
	 vIt != markstats->GetValidLabelValues().end();
	 ++vIt){
	if ( markstats->HasLabel(*vIt) ){
	    typename OutputImType::PixelType labelValue = *vIt;
	    if ( markstats->GetCount(labelValue) > msize){
		msize = markstats->GetCount(labelValue);
		MMean = markstats->GetMean(labelValue);
		MSigma = markstats->GetSigma(labelValue);
		}
	    }
	}

    std::cout << "Stats " << msize << " " << MMean << " " << MSigma << std::endl;

    // create a mask to remove bright border pixels
    typedef typename itk::Image<unsigned char, InImType::ImageDimension > MaskImType;
    typedef typename itk::BinaryThresholdImageFilter<InImType, MaskImType> ThreshType;
    typename ThreshType::Pointer thresh = ThreshType::New();
    thresh->SetInput(input);
    thresh->SetUpperThreshold(MMean - MSigma);
    thresh->SetInsideValue(1);
    thresh->SetOutsideValue(0);

    typedef typename itk::MaskImageFilter<OutputImType, MaskImType, OutputImType> MaskerType;
    typename MaskerType::Pointer masker = MaskerType::New();
    masker->SetInput(select->GetOutput());
    masker->SetInput2(thresh->GetOutput());

    typename OutputImType::Pointer result = masker->GetOutput();
    result->Update();
    result->DisconnectPipeline();
    return(result);
    }


template <class InImType, class OutputImType> 
typename OutputImType::Pointer findDarkMarkers(
    typename InImType::Pointer input,
    int radius,
    int offset){

  typename InImType::Pointer bthI = 0;
  {
    // Now for traditional dark markers 
  typedef typename itk::FlatStructuringElement< InImType::ImageDimension > SRType;
  typename SRType::RadiusType rad, bthrad;
  rad.Fill(radius);
  SRType kernel, bthkernel;
  
  kernel = SRType::Box(rad);
  
  typedef typename itk::GrayscaleDilateImageFilter<InImType, InImType, SRType> DilateType;
  typename DilateType::Pointer dilate = DilateType::New();
  dilate->SetInput(input);
  dilate->SetKernel(kernel);

  // black tophat filter - there is some sort of edge artifact -
  // uneven brightness
  typedef typename itk::BlackTopHatImageFilter<InImType, InImType, SRType> BTHType;
  typename BTHType::Pointer bth = BTHType::New();
  bthrad.Fill(radius*10);
  bthkernel = SRType::Box(bthrad);
  bth->SetInput(dilate->GetOutput());
  bth->SetKernel(bthkernel);
  bthI = bth->GetOutput();
  bthI->Update();
  bthI->DisconnectPipeline();
  }
  typedef typename itk::Image<unsigned char, InImType::ImageDimension > BinaryType;
  typedef typename itk::OtsuThresholdImageFilter<InImType, OutputImType> OtsuType;
  typename OtsuType::Pointer Otsu = OtsuType::New();
  Otsu->SetInput(bthI);
  Otsu->SetInsideValue(0);
  Otsu->SetOutsideValue(offset);

  // not sure if we need to worry about markerssofar - could use it to
  // delete some. Instead we'll take the max later
  typedef typename itk::ConnectedComponentImageFilter<OutputImType, OutputImType> LabellerType;
  typename LabellerType::Pointer labeller = LabellerType::New();
  labeller->SetInput(Otsu->GetOutput());
  labeller->SetBackgroundValue(0);

  typedef typename itk::AddImageFilter<OutputImType, OutputImType, OutputImType> AdderType;
  typename AdderType::Pointer adder = AdderType::New();
  adder->SetInput(Otsu->GetOutput());
  adder->SetInput2(labeller->GetOutput());
  
  
  typename OutputImType::Pointer result = adder->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  writeImDbg<OutputImType>(result, "mk.mha");
  return(result);

}

template <class InImType> 
typename InImType::Pointer computeGrad(
    typename InImType::Pointer input,
    int radius,
    float smooth){
    // Compute a smoothed morphological gradient - original - erosion
    // optionally smooth it. This gradient will put the peaks on the
    // outside of a thin dark line, which I think is desirable for this
    // application.
    typedef typename itk::FlatStructuringElement< InImType::ImageDimension > SRType;
    typename SRType::RadiusType rad;
    rad.Fill(radius);
    SRType kernel;

    kernel = SRType::Box(rad);

    typedef typename itk::GrayscaleErodeImageFilter<InImType, InImType, SRType> ErodeType;
    typename ErodeType::Pointer erode = ErodeType::New();
    erode->SetInput(input);
    erode->SetKernel(kernel);

    typedef typename itk::SubtractImageFilter<InImType, InImType, InImType> SubType;
    typename SubType::Pointer sub = SubType::New();
    sub->SetInput(input);
    sub->SetInput2(erode->GetOutput());
    typename InImType::Pointer result = sub->GetOutput();
    typedef typename itk::SmoothingRecursiveGaussianImageFilter<InImType, InImType> SmoothType;
    typename SmoothType::Pointer smoother = SmoothType::New();
  
    typename InImType::SpacingType sp = input->GetSpacing();

    if (smooth > 0){
	// Lazy - assume isotropic
	float sigma = sp[0] * smooth;
	smoother->SetInput(sub->GetOutput());
	smoother->SetSigma(sigma);
	result = smoother->GetOutput();
	}

    result->Update();
    result->DisconnectPipeline();
    return(result);
    }
#define USE_FLOAT
#define USE_UI32
template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(CmdLineType &CmdLineObj){
#ifdef USE_FLOAT
    typedef float  TRealType;
    std::cerr << "Using single precision (float)." << std::endl;
#else
    typedef double TRealType;
    std::cerr << "Using double precision (double)." << std::endl;
#endif

#ifdef USE_UI32
    typedef uint32_t  OutputPixelType;
    std::cerr << "Using uInt32." << std::endl;
#else
    typedef uint64_t  OutputPixelType;
    std::cerr << "Using uInt64." << std::endl;
#endif
    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<TRealType, Dimension>        GreyImageType;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;

    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(CmdLineObj.InputImFile);
    //reader->ReleaseDataFlagOn();
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

    
    std::cout << "Bright markers" << std::endl;
    // Markers will be large areas
    int ForegroundLabels = 0;
    typename OutputImageType::Pointer brightmarkers =
      findBrightMarkers<InputImageType, OutputImageType>(reader->GetOutput(), CmdLineObj.markererode, CmdLineObj.markervolume, ForegroundLabels);

    std::cout << "Dark markers" << std::endl;
    typename OutputImageType::Pointer darkmarkers = findDarkMarkersWS<InputImageType, OutputImageType>(reader->GetOutput(), brightmarkers, ForegroundLabels+1);

    std::cout << "Dark markers blood cells" << std::endl;

    typename OutputImageType::Pointer darkmarkersBC = findDarkMarkers<InputImageType, OutputImageType>(reader->GetOutput(), CmdLineObj.darkmarkerdilate, ForegroundLabels+2);

    typedef typename itk::MaximumImageFilter<OutputImageType, OutputImageType, OutputImageType> MaxType;
    typename MaxType::Pointer combdark = MaxType::New();
    combdark->SetInput(darkmarkers);
    combdark->SetInput2(darkmarkersBC);

    typename MaxType::Pointer comb = MaxType::New();
    comb->SetInput(brightmarkers);
    comb->SetInput2(combdark->GetOutput());
       
    typename OutputImageType::Pointer finalmarkers = comb->GetOutput();
    finalmarkers->Update();
    finalmarkers->DisconnectPipeline();
    // dispose intermediate steps
    brightmarkers = 0;
    darkmarkers = 0;
    darkmarkersBC = 0;
    comb = 0;
    combdark = 0;

    std::cout << "Gradient" << std::endl;
    typename InputImageType::Pointer grad = computeGrad<InputImageType>(reader->GetOutput(), CmdLineObj.MorphGradRad, CmdLineObj.SmoothGradSigma);

    reader = 0;
    writeImDbg<OutputImageType>(finalmarkers, "markers.mha");

    std::cout << "Watershed" << std::endl;

    typedef typename itk::MorphologicalWatershedFromMarkersImageFilter<InputImageType, OutputImageType> WSType;
    typename WSType::Pointer ws = WSType::New();
    FilterWatcher watcherWS2(ws);
    ws->SetInput(grad);
    ws->SetMarkerImage(finalmarkers);
    ws->SetMarkWatershedLine(false);


    typedef typename itk::ChangeLabelImageFilter<OutputImageType, OutputImageType> ChangeType;
    typename ChangeType::Pointer changer = ChangeType::New();
    changer->SetInput(ws->GetOutput());
    changer->SetChange(ForegroundLabels+1, 0);

    writeIm<OutputImageType>(changer->GetOutput(), CmdLineObj.OutputImFile);


    return EXIT_SUCCESS;
 

}


template<typename InputComponentType, typename InputPixelType>
int dispatch_D(size_t dimensionType, CmdLineType &CmdLineObj){
    int res= 0;
    switch (dimensionType){
    // case 1:
    //     res= DoIt<InputComponentType, InputPixelType, 1>(argc, argv);
    //     break;
    case 2:
        res= DoIt<InputComponentType, InputPixelType, 2>(CmdLineObj);
        break;
    case 3:
        res= DoIt<InputComponentType, InputPixelType, 3>(CmdLineObj);
        break;
    default:
        std::cerr << "Error: Images of dimension " << dimensionType << " are not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

template<typename InputComponentType>
int dispatch_pT(itk::ImageIOBase::IOPixelType pixelType, size_t dimensionType, CmdLineType &CmdLineObj){
    int res= 0;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (pixelType){
    case itk::ImageIOBase::SCALAR:{
        typedef InputComponentType InputPixelType;
        res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, CmdLineObj);
        } break;
    case itk::ImageIOBase::UNKNOWNPIXELTYPE:
    default:
        std::cerr << std::endl << "Error: Pixel type not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

int dispatch_cT(itk::ImageIOBase::IOComponentType componentType, itk::ImageIOBase::IOPixelType pixelType, size_t dimensionType, CmdLineType &CmdLineObj){
    int res= 0;

    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#a8dc783055a0af6f0a5a26cb080feb178
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00107
    //IOComponentType: UNKNOWNCOMPONENTTYPE, UCHAR, CHAR, USHORT, SHORT, UINT, INT, ULONG, LONG, FLOAT, DOUBLE

    switch (componentType){
    case itk::ImageIOBase::UCHAR:{        // uint8_t
        typedef unsigned char InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, CmdLineObj);
        } break;
    case itk::ImageIOBase::CHAR:{         // int8_t
        typedef char InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, CmdLineObj);
        } break;
    case itk::ImageIOBase::USHORT:{       // uint16_t
        typedef unsigned short InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, CmdLineObj);
        } break;
    case itk::ImageIOBase::SHORT:{        // int16_t
        typedef short InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, CmdLineObj);
        } break;
    case itk::ImageIOBase::UINT:{         // uint32_t
        typedef unsigned int InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, CmdLineObj);
        } break;
    case itk::ImageIOBase::INT:{          // int32_t
        typedef int InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, CmdLineObj);
        } break;
    case itk::ImageIOBase::ULONG:{        // uint64_t
        typedef unsigned long InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, CmdLineObj);
        } break;
    case itk::ImageIOBase::LONG:{         // int64_t
        typedef long InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, CmdLineObj);
        } break;
    case itk::ImageIOBase::FLOAT:{        // float32
        typedef float InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, CmdLineObj);
        } break;
    case itk::ImageIOBase::DOUBLE:{       // float64
        typedef double InputComponentType;
        res= dispatch_pT<InputComponentType>(pixelType, dimensionType, CmdLineObj);
        } break;
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
    default:
        std::cerr << "unknown component type" << std::endl;
        break;
        }//switch
    return res;
    }


void GetImageType (std::string fileName,
    itk::ImageIOBase::IOPixelType &pixelType,
    itk::ImageIOBase::IOComponentType &componentType,
    size_t &dimensionType
    ){
    typedef itk::Image<char, 1> ImageType; //template initialization parameters need to be given but can be arbitrary here
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

  // tclap processing - see if you like it.

    itk::ImageIOBase::IOPixelType pixelType;
    typename itk::ImageIOBase::IOComponentType componentType;
    size_t dimensionType;

    CmdLineType CmdLineObj;
    ParseCmdLine(argc, argv, CmdLineObj);


    try {
    GetImageType(CmdLineObj.InputImFile, pixelType, componentType, dimensionType);
        }//try
    catch( itk::ExceptionObject &excep){
        std::cerr << argv[0] << ": exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
        }

    return dispatch_cT(componentType, pixelType, dimensionType,CmdLineObj);
    }

