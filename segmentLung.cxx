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
#include <itkThresholdImageFilter.h>
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
#include <itkChangeLabelLabelMapFilter.h>
#include <itkAddImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkAttributeOpeningLabelMapFilter.h>
#include <itkRegionalMaximaImageFilter.h>
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
  int DarkVol, markererode, markervolume, darkmarkervolume, MorphGradRad, darkmarkerdilate, BGCorrectionRadius;
  float SmoothGradSigma, darksmooth;
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

//	ValueArg<int> darkmarkerdilateArg("","darkmarkerdilate","Radius of dilation of image used to produce bloodcell markers (voxels)",true, 20, "integer");
//	cmd.add( darkmarkerdilateArg );
	ValueArg<int> bgcorrectArg("","bgcorrect","Radius of tophat filter used in background illumination correction",true, 50, "integer");
	cmd.add( bgcorrectArg );

	ValueArg<int> markervolArg("","markervol","Minimum volume of a bright marker",false, 100000, "integer");
	cmd.add( markervolArg);

//	ValueArg<int> darkmarkervolArg("","darkmarkervol","Maximum volume of a dark marker",false, 1000, "integer");
//	cmd.add( darkmarkervolArg);

	ValueArg<int> morphgradArg("","morphgradrad","Radius of SE for morphologicad gradient",false, 1, "integer");
	cmd.add( morphgradArg );

	ValueArg<float> gradsmoothArg("","gradsigma","Smoothing kernel (voxels) in gradient calculation",false, 1, "float");
	cmd.add( gradsmoothArg);

	ValueArg<float> darksmoothArg("","darksigma","Smoothing kernel (voxels) for regional minima (dark markers)",false, 5, "float");
	cmd.add( darksmoothArg);

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
	CmdLineObj.BGCorrectionRadius = bgcorrectArg.getValue();
	//CmdLineObj.darkmarkervolume = darkmarkervolArg.getValue();
	CmdLineObj.MorphGradRad = morphgradArg.getValue();
	CmdLineObj.SmoothGradSigma = gradsmoothArg.getValue();
	CmdLineObj.darksmooth = darksmoothArg.getValue();
	//CmdLineObj.darkmarkerdilate = darkmarkerdilateArg.getValue();
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
    Otsu->SetInsideValue(1);
    Otsu->SetOutsideValue(0);

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

template <class InImType, class OutputImType, class BinaryImType> 
typename BinaryImType::Pointer findDarkMarkersWS(
    typename InImType::Pointer input,
    typename OutputImType::Pointer markers,
    int BorderLabel){
   
  // assumes that input is a blacktophat filter output
    // typedef typename itk::StatisticsImageFilter<InImType> StatsType;
    // typename StatsType::Pointer stats = StatsType::New();
    // stats->SetInput(input);
    // stats->Update();

    // typedef typename itk::ShiftScaleImageFilter<InImType, InImType> ShiftType;
    // typename ShiftType::Pointer shifter = ShiftType::New();
    // shifter->SetInput(input);
    // shifter->SetShift(- stats->GetMaximum() );
    // shifter->SetScale(-1);

    typedef typename itk::MorphologicalWatershedFromMarkersImageFilter<InImType, OutputImType> WSType;
    typename WSType::Pointer ws = WSType::New();
    FilterWatcher watcherws1(ws);
    ws->SetInput(input);
    ws->SetMarkerImage(markers);
    ws->SetMarkWatershedLine(true);

    typedef typename itk::BinaryThresholdImageFilter<OutputImType, BinaryImType> SelectType;
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

    // create a mask to remove dark border pixels
    //typedef typename itk::Image<unsigned char,
    //InImType::ImageDimension > MaskImType;
    
    typedef typename itk::BinaryThresholdImageFilter<InImType, BinaryImType> ThreshType;
    typename ThreshType::Pointer thresh = ThreshType::New();
    thresh->SetInput(input);
    thresh->SetUpperThreshold(MMean + MSigma);
    thresh->SetInsideValue(0);
    thresh->SetOutsideValue(1);

    typedef typename itk::MaskImageFilter<BinaryImType, BinaryImType, BinaryImType> MaskerType;
    typename MaskerType::Pointer masker = MaskerType::New();
    masker->SetInput(select->GetOutput());
    masker->SetInput2(thresh->GetOutput());

    typename BinaryImType::Pointer result = masker->GetOutput();
    result->Update();
    result->DisconnectPipeline();
    return(result);
    }

///////////////////////////////////////////////////////////////////////////////
template <class InImType>
typename InImType::Pointer adjustBG(
    typename InImType::Pointer input,
    int BTHradius) {
// BTHradius - size of filter used to correct background (may not be

  typedef typename itk::FlatStructuringElement< InImType::ImageDimension > SRType;
  typename SRType::RadiusType rad, bthrad;
  rad.Fill(BTHradius);
  SRType kernel, bthkernel;
  
  kernel = SRType::Box(rad);

  // black tophat filter - there is some sort of edge artifact -
  // uneven brightness
  typedef typename itk::BlackTopHatImageFilter<InImType, InImType, SRType> BTHType;
  typename BTHType::Pointer bth = BTHType::New();
  bthrad.Fill(BTHradius);
  bthkernel = SRType::Box(bthrad);
  bth->SetInput(input);
  bth->SetKernel(bthkernel);
  typename InImType::Pointer bthI = bth->GetOutput();
  bthI->Update();
  bthI->DisconnectPipeline();
  return(bthI);
}

///////////////////////////////////////////////////////////////////////////////
template <class InImType, class OutputImType> 
typename OutputImType::Pointer findDarkMarkersMin(
    typename InImType::Pointer input,
    float smooth,
    int offset) {

// Take 2 of finding dark markers. Some red blood cells are very close
// to the other dark structures, and aren't separated by simple
// procedures, like those used in Take 1 (below). Thus we'll try
// something more like a regional minima approach, but with some
// tweaks to avoid the need for preflooding, and tricks to make it
// faster.
// Still need to adjust background

// necessary in this approach.
// smooth - smoothing kernel.
// offset - added to the label values
  typedef typename itk::Image<unsigned char, InImType::ImageDimension > BinaryType;

  typename InImType::Pointer bthI = 0;
  typename BinaryType::Pointer darkminI = 0;
  {

  typedef typename itk::SmoothingRecursiveGaussianImageFilter<InImType, InImType> SmoothType;
  typename SmoothType::Pointer smoother = SmoothType::New();
  
  typename InImType::SpacingType sp = input->GetSpacing();
  
  if (smooth > 0){
  // Lazy - assume isotropic
  float sigma = sp[0] * smooth;
  smoother->SetInput(input);
  smoother->SetSigma(sigma);
  bthI = smoother->GetOutput();
  } else {
  bthI = input;
  }
  bthI->Update();
  bthI->DisconnectPipeline();
  }
  // Now we want to compute a threshold, and only look for dark makers
  // in the dark part. However the black top hat has swapped things
  // around
  writeImDbg<InImType>(bthI, "smoothedBTH.mha");
  {
  // Use the otsu threshold filter, even though we only want the
  // threshold value.
  typedef typename itk::OtsuThresholdImageFilter<InImType, BinaryType> OtsuType;
  typename OtsuType::Pointer Otsu = OtsuType::New();
  Otsu->SetInput(bthI);
  Otsu->SetInsideValue(0);
  Otsu->SetOutsideValue(1);
  Otsu->Update();

  typename InImType::PixelType othresh = Otsu->GetThreshold();

  // This is to set the areas we aren't interested to a value that the
  // regional extrema filter will ignore.
  typedef typename itk::ThresholdImageFilter<InImType> ThresholdType;
  typename ThresholdType::Pointer tt = ThresholdType::New();
  tt->SetInput(bthI);
  tt->SetOutsideValue(itk::NumericTraits< typename InImType::PixelType >::NonpositiveMin() );
  tt->ThresholdBelow(othresh);
  
  typedef typename itk::RegionalMaximaImageFilter<InImType, BinaryType> RegMaxType;
  typename RegMaxType::Pointer rmax = RegMaxType::New();
  FilterWatcher watcherrmax(rmax);
  rmax->SetForegroundValue(1);
  rmax->SetBackgroundValue(0);
  rmax->SetInput(tt->GetOutput());
  darkminI = rmax->GetOutput();
  darkminI->Update();
  darkminI->DisconnectPipeline();
  writeImDbg<BinaryType>(darkminI, "regmin.mha");

  }
  return(darkminI);
#if 0
  // label maps to reduce memory footprint
  typedef typename itk::BinaryImageToShapeLabelMapFilter<BinaryType> LabellerType;
  typename LabellerType::Pointer labeller = LabellerType::New();
  labeller->SetInput(darkminI);
  labeller->SetInputForegroundValue(1);
  // add the offset
  typename LabellerType::OutputImageType::Pointer rleObj = labeller->GetOutput();
  rleObj->Update();
  rleObj->DisconnectPipeline();

  for(unsigned int i = 1; i < rleObj->GetNumberOfLabelObjects(); ++i)
   {
   typename LabellerType::OutputImageType::LabelObjectType* shapeLabelObject =
     rleObj->GetNthLabelObject(i);
     shapeLabelObject->SetLabel(i + offset + 1);
   }

  typedef typename itk::LabelMapToLabelImageFilter<typename LabellerType::OutputImageType, OutputImType> ConvType;

  typename ConvType::Pointer tolabIm = ConvType::New();
  tolabIm->SetInput(rleObj);
  
  typename OutputImType::Pointer result = tolabIm->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  writeImDbg<OutputImType>(result, "mk.mha");
  return(result);
#endif

}
//////////////////////////////////////////////////////////////////////////////
template <class InImType, class OutputImType> 
typename OutputImType::Pointer findDarkMarkers(
    typename InImType::Pointer input,
    int radius,
    int offset,
    float volumeopening){

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
  typedef typename itk::OtsuThresholdImageFilter<InImType, BinaryType> OtsuType;
  typename OtsuType::Pointer Otsu = OtsuType::New();
  Otsu->SetInput(bthI);
  Otsu->SetInsideValue(0);
  Otsu->SetOutsideValue(1);

  // retain small regions
  typedef typename itk::BinaryImageToShapeLabelMapFilter<BinaryType> LabellerType;
  typename LabellerType::Pointer labeller = LabellerType::New();
  labeller->SetInput(Otsu->GetOutput());
  labeller->SetInputForegroundValue(1);
  typedef typename itk::ShapeOpeningLabelMapFilter<typename LabellerType::OutputImageType> ShapeFilterType;
  typename ShapeFilterType::Pointer shapefilter = ShapeFilterType::New();
  shapefilter->SetLambda(volumeopening);
  shapefilter->SetAttribute("NumberOfPixels");
  shapefilter->ReverseOrderingOn();
  shapefilter->SetInput(labeller->GetOutput());

  typename LabellerType::OutputImageType::Pointer rleObj = shapefilter->GetOutput();
  rleObj->Update();
  rleObj->DisconnectPipeline();

  for(unsigned int i = 1; i < rleObj->GetNumberOfLabelObjects(); ++i)
   {
   typename LabellerType::OutputImageType::LabelObjectType* shapeLabelObject =
     rleObj->GetNthLabelObject(i);
     shapeLabelObject->SetLabel(i + offset + 1);
   }
 


  typedef typename itk::LabelMapToLabelImageFilter<typename LabellerType::OutputImageType, OutputImType> ConvType;

  typename ConvType::Pointer tolabIm = ConvType::New();
  tolabIm->SetInput(rleObj);
  
  typename OutputImType::Pointer result = tolabIm->GetOutput();
  result->Update();
  result->DisconnectPipeline();
  writeImDbg<OutputImType>(result, "mk.mha");
  return(result);


#if 0
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
#endif
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
    typedef typename itk::Image<unsigned char, Dimension > BinaryImType;

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

    std::cout << "Background subtraction" << std::endl;
    typename InputImageType::Pointer bgcorrected = adjustBG<InputImageType>(reader->GetOutput(), CmdLineObj.BGCorrectionRadius);
    
    writeImDbg<InputImageType>(bgcorrected, "bth.mha");

    std::cout << "Bright markers" << std::endl;
    // Markers will be large areas
    int ForegroundLabels = 0;
    typename OutputImageType::Pointer brightmarkers =
      findBrightMarkers<InputImageType, OutputImageType>(bgcorrected, CmdLineObj.markererode, CmdLineObj.markervolume, ForegroundLabels);

    std::cout << "Dark markers" << std::endl;
    typename BinaryImType::Pointer darkmarkers = findDarkMarkersWS<InputImageType, OutputImageType, BinaryImType>(bgcorrected, brightmarkers, 1);

    std::cout << "Dark markers blood cells" << std::endl;

    typename BinaryImType::Pointer darkmarkersBC = findDarkMarkersMin<InputImageType, BinaryImType>(bgcorrected, CmdLineObj.darksmooth, ForegroundLabels+2);
    typename OutputImageType::Pointer alldarkmarkers = 0; 

    std::cout << "Dark marker combination" << std::endl;

    typedef typename itk::MaximumImageFilter<OutputImageType, OutputImageType, OutputImageType> MaxType;
    typedef typename itk::MaximumImageFilter<BinaryImType, BinaryImType, BinaryImType> BMaxType;
    {
    typename BMaxType::Pointer combdark = BMaxType::New();
    combdark->SetInput(darkmarkers);
    combdark->SetInput2(darkmarkersBC);

    // now we label the combined markers
    typedef typename itk::ConnectedComponentImageFilter<BinaryImType, OutputImageType> LabellerType;
    typename LabellerType::Pointer labeller = LabellerType::New();
    labeller->SetInput(combdark->GetOutput());
    labeller->SetFullyConnected(true);
    labeller->SetBackgroundValue(0);

    // add the offset and mask
    typedef typename itk::AddImageFilter<OutputImageType, OutputImageType, OutputImageType> AddType;
    typename AddType::Pointer adder = AddType::New();
    adder->SetInput1(labeller->GetOutput());
    adder->SetConstant2(ForegroundLabels + 1);

    typedef typename itk::MaskImageFilter<OutputImageType, BinaryImType, OutputImageType> MaskerType;
    typename MaskerType::Pointer masker = MaskerType::New();
    masker->SetInput(adder->GetOutput());
    masker->SetInput2(combdark->GetOutput());
    alldarkmarkers = masker->GetOutput();
    alldarkmarkers->Update();
    alldarkmarkers->DisconnectPipeline();
    std::cout << "Alldarkmarkers done" << std::endl;
    }

    typename MaxType::Pointer comb = MaxType::New();
    comb->SetInput(brightmarkers);
    comb->SetInput2(alldarkmarkers);
       
    typename OutputImageType::Pointer finalmarkers = comb->GetOutput();
    finalmarkers->Update();
    finalmarkers->DisconnectPipeline();
    // dispose intermediate steps
    brightmarkers = 0;
    darkmarkers = 0;
    darkmarkersBC = 0;
    comb = 0;
    //combdark = 0;

    std::cout << "Gradient" << std::endl;
    typename InputImageType::Pointer grad = computeGrad<InputImageType>(reader->GetOutput(), CmdLineObj.MorphGradRad, CmdLineObj.SmoothGradSigma);
    writeImDbg<InputImageType>(grad, "grad.mha");
    reader = 0;
    writeImDbg<OutputImageType>(finalmarkers, "markers.mha");

    std::cout << "Watershed" << std::endl;

    typedef typename itk::MorphologicalWatershedFromMarkersImageFilter<InputImageType, OutputImageType> WSType;
    typename WSType::Pointer ws = WSType::New();
    FilterWatcher watcherWS2(ws);
    ws->SetInput(grad);
    ws->SetMarkerImage(finalmarkers);
    ws->SetMarkWatershedLine(false);

    // typedef typename itk::ChangeLabelImageFilter<OutputImageType, OutputImageType> ChangeType;
    // writeImDbg<OutputImageType>(ws->GetOutput(), "wsout.mha");

    // typename ChangeType::Pointer changer = ChangeType::New();
    // changer->SetInput(ws->GetOutput());
    // changer->SetChange(ForegroundLabels+1, 0);

    writeIm<OutputImageType>(ws->GetOutput(), CmdLineObj.OutputImFile);

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

