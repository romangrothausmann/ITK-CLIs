#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkShiftScaleImageFilter.h>
#include <itkMorphologicalWatershedImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkAreaClosingImageFilter.h>

#include "tclap/CmdLine.h"

bool WriteDebug = false;

template <class TImage>
void writeIm(typename TImage::Pointer Im, std::string filename)
{
  typedef typename itk::ImageFileWriter<TImage> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput(Im);
  writer->SetFileName(filename.c_str());
  writer->Update();
}


template <class ImType>
void writeImDbg(typename ImType::Pointer Im, std::string filename)
{
  if (WriteDebug)
    {
    writeIm<ImType>(Im, filename);
    }
}


typedef class CmdLineType
{
public:
  std::string InputImFile, OutputImFile;
  int DarkVol;
  bool Compression;
} CmdLineType;

void ParseCmdLine(int argc, char* argv[],
                  CmdLineType &CmdLineObj
                  )
{
  using namespace TCLAP;
  try
    {
    // Define the command line object.
    CmdLine cmd("Lung segmentation prototype ", ' ', "0.9");


    ValueArg<std::string> inArg("i","input","input image",true,"result","string");
    cmd.add( inArg );

    ValueArg<std::string> outArg("o","output","output image", true,"","string");
    cmd.add( outArg );

    ValueArg<int> darkArg("","darkvol","Volume used for attribute closing (voxels)",true, 3000, "integer");
    cmd.add( darkArg);

    SwitchArg compArg("c", "compress", "write compressed images", false);
    cmd.add( compArg );

   // Parse the args.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg.
    CmdLineObj.InputImFile = inArg.getValue();
    CmdLineObj.OutputImFile = outArg.getValue();
    CmdLineObj.DarkVol = darkArg.getValue();
    CmdLineObj.Compression = compArg.getValue();

    }
  catch (ArgException &e)  // catch any exceptions
    {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

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
    reader->ReleaseDataFlagOn();
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

    // Try some preprocessing - attribute filter
    typedef typename itk::AreaClosingImageFilter<InputImageType, InputImageType> ACType;
    typename ACType::Pointer areaclose = ACType::New();
    areaclose->SetInput(reader->GetOutput());
    areaclose->SetLambda(CmdLineObj.DarkVol);
    areaclose->UseImageSpacingOff();

    typedef itk::ImageFileWriter<InputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    
    FilterWatcher watcherO(writer);
    writer->SetFileName(CmdLineObj.OutputImFile);
    writer->SetInput(areaclose->GetOutput());
    writer->SetUseCompression(CmdLineObj.Compression);
    try{
        writer->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

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
        GetImageType(argv[1], pixelType, componentType, dimensionType);
        }//try
    catch( itk::ExceptionObject &excep){
        std::cerr << argv[0] << ": exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
        }

    return dispatch_cT(componentType, pixelType, dimensionType,CmdLineObj);
    }

