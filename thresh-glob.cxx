////program for itkBinaryThresholdImageFilter
//01: based on template.cxx


#include <complex>

#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkBinaryThresholdImageFilter.h>
#ifdef USE_SDI
#include <itkPipelineMonitorImageFilter.h>
#endif


int dispatch_cT(itk::ImageIOBase::IOPixelType, itk::ImageIOBase::IOComponentType, size_t, int, char **);

template<typename InputComponentType>
int dispatch_pT(itk::ImageIOBase::IOPixelType pixelType, size_t, int, char **);

template<typename InputComponentType, typename InputPixelType>
int dispatch_D(size_t, int, char **);

template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int, char *argv[]);






template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){

    typedef InputPixelType  OutputPixelType; #needs more sophisticated check

    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;


    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(argv[1]);
#ifndef USE_SDI
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
#endif

    typename InputImageType::Pointer input= reader->GetOutput();



    typedef itk::BinaryThresholdImageFilter<InputImageType, OutputImageType> FilterType;
    typename FilterType::Pointer filter= FilterType::New();
    filter->SetInput(input);

    InputPixelType th_l, th_h; //not InputPixelType!


    char* arg;

    arg= argv[4];
    if (!strcasecmp(arg, "max")){
        fprintf(stderr, "%s: identified as maximum\n", arg);
        th_l= itk::NumericTraits<InputPixelType>::max();
        //itk::NumericTraits<InputPixelType>::max(&pv);
        std::cerr << "Max: " << +th_l << std::endl;
        }
    else if (!strcasecmp(arg, "min")){
        fprintf(stderr, "%s: identified as minimum\n", arg);
        th_l= itk::NumericTraits<InputPixelType>::min();
        std::cerr << "Min: " << +th_l << std::endl;
        std::cerr << "Is -1 negative? " << itk::NumericTraits<InputPixelType>::IsNegative(-1) << std::endl;
        std::cerr << "Is 1 negative? " << itk::NumericTraits<InputPixelType>::IsNegative(1) << std::endl;
        }
    else
        th_l= InputPixelType(atof(arg));

    arg= argv[5];
    if (!strcasecmp(arg, "max")){
        fprintf(stderr, "%s: identified as maximum\n", arg);
        th_h= itk::NumericTraits<InputPixelType>::max();
        }
    else if (!strcasecmp(arg, "min")){
        fprintf(stderr, "%s: identified as minimum\n", arg);
        th_h= itk::NumericTraits<InputPixelType>::min();
        }
    else
        th_h= InputPixelType(atof(arg));

    if (th_l > th_h){
        filter->SetLowerThreshold(th_h);
        filter->SetUpperThreshold(th_l);
	if(argc > 6)
	    filter->SetOutsideValue(static_cast<OutputPixelType>(atof(argv[6])));
	else
	    filter->SetOutsideValue(itk::NumericTraits<OutputPixelType>::max());
        filter->SetInsideValue(itk::NumericTraits<OutputPixelType>::Zero);
        }
    else {
        filter->SetLowerThreshold(th_l);
        filter->SetUpperThreshold(th_h);
        filter->SetOutsideValue(itk::NumericTraits<OutputPixelType>::Zero);
	if(argc > 6)
	    filter->SetInsideValue(static_cast<OutputPixelType>(atof(argv[6])));
	else
	    filter->SetInsideValue(itk::NumericTraits<OutputPixelType>::max());
        }

    std::cerr << "lower_th: "<< +filter->GetLowerThreshold() << "   upper_th: " << +filter->GetUpperThreshold() << std::endl; //+ promotes variable to a type printable as a number (e.g. for char)

#ifndef USE_SDI
    FilterWatcher watcher1(filter);
    try{
        filter->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }
#else
    typedef itk::PipelineMonitorImageFilter<OutputImageType> MonitorFilterType;
    typename MonitorFilterType::Pointer monitorFilter = MonitorFilterType::New();
    monitorFilter->SetInput(filter->GetOutput());
#endif


    typename OutputImageType::Pointer output= filter->GetOutput();

    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer);
    writer->SetFileName(argv[2]);
#ifndef USE_SDI
    writer->SetInput(output);
    writer->SetUseCompression(atoi(argv[3]));
#else
    writer->SetInput(monitorFilter->GetOutput());
    writer->UseCompressionOff(); //writing compressed is not supported for streaming!
    writer->SetNumberOfStreamDivisions(atoi(argv[3]));
#endif
    try{
        writer->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

#ifdef USE_SDI
    if (!monitorFilter->VerifyAllInputCanStream(atoi(argv[3]))){
        //std::cerr << monitorFilter;
        }
#endif

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
    case 1:
        res= DoIt<InputComponentType, InputPixelType, 1>(argc, argv);
        break;
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
    if ( argc < 6 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Output_Image"
#ifndef USE_SDI
                  << " compress"
#else
                  << " stream-chunks"
#endif
		  << " lower upper"
		  << " [inside-value]"
		  << std::endl;

        return EXIT_FAILURE;
        }

    if(atoi(argv[3]) < 0){
        std::cerr << "3rd parameter must not be negative" << std::endl;
        return EXIT_FAILURE;
        }

#ifndef USE_SDI
    if(atoi(argv[3]) > 1){
        std::cerr << "compress must be 0 or 1 (to avoid confusion with stream-chunks of SDI version" << std::endl;
        return EXIT_FAILURE;
        }
#else
    if(atoi(argv[3]) < 2){
        std::cerr << "stream-chunks must be > 1 (or use non-SDI version)" << std::endl;
        return EXIT_FAILURE;
        }
#endif


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






