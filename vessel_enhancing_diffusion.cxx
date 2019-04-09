////program for itkVesselEnhancingDiffusion3DImageFilter
///? a rewrite of (Enquobahrie2007) Enquobahrie, A.; Ibáñez, L.; Bullitt, E.; Aylward, S. (2007), which was never included in ITK? but ported to ITK4: https://github.com/rcasero/gerardus/tree/master/cpp/src/third-party/IJ-Vessel_Enhancement_Diffusion.1
///? based on (Manniesing2006) Manniesing, R.; Viergever, M. A.; Niessen, W. J. (2006) Vessel enhancing diffusion: A scale space representation of vessel structures
//01: based on template.cxx, vtkVmtk/Segmentation/vtkvmtkVesselEnhancingDiffusion3DImageFilter.h and https://github.com/InsightSoftwareConsortium/LesionSizingToolkit/blob/master/test/itkVEDTest.cxx


#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkVesselEnhancingDiffusion3DImageFilter.h>//ITK-4.7.0: Module LesionSizingToolkit
#include <itkImageFileWriter.h>

/*
enum{
    EQUISPACED_STEPS,
    LOGARITHMIC_STEPS
    };

double ComputeSigmaValue(int scaleLevel, int SigmaStepMethod){
    double sigmaValue;

    if (this->NumberOfSigmaSteps < 2)
	return this->SigmaMin;

    switch (SigmaStepMethod){
    case EQUISPACED_STEPS:{
        double stepSize = (SigmaMax - SigmaMin) / (NumberOfSigmaSteps - 1);
        if (stepSize < 1e-10)
	    stepSize = 1e-10;
        sigmaValue = SigmaMin + stepSize * scaleLevel;
	} break;
    case LOGARITHMIC_STEPS:{
	double stepSize = (vcl_log(SigmaMax) - vcl_log(SigmaMin)) / (NumberOfSigmaSteps - 1);
	if (stepSize < 1e-10)
	    stepSize = 1e-10;
	sigmaValue = vcl_exp(vcl_log (SigmaMin) + stepSize * scaleLevel);
	} break;
    default:
	vtkErrorMacro("Error: undefined sigma step method.");
	sigmaValue = 0.0;
	break;
	}

    return sigmaValue;
    }
*/

template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){

    typedef InputPixelType  OutputPixelType;

    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;

    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(argv[1]);
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

    const typename InputImageType::Pointer& input= reader->GetOutput();

    const double SigmaMin= atof(argv[4]);
    const double SigmaMax= atof(argv[5]);
    const int NumberOfSigmaSteps= atoi(argv[6]);

    if (NumberOfSigmaSteps < 2){
	std::cerr << "sigmaSteps has to be greater than 1!" << std::endl;
	return EXIT_FAILURE;
	}


    typedef itk::VesselEnhancingDiffusion3DImageFilter<InputPixelType, Dimension> FilterType;
    typename FilterType::Pointer filter= FilterType::New();
    filter->SetInput(input);
    filter->ReleaseDataFlagOn();
    filter->SetDefaultPars();//sets default values
    // filter->SetAlpha();
    // filter->SetBeta();
    // filter->SetGamma();
    // filter->SetEpsilon();
    // filter->SetOmega();//? AnisotropicDiffusionVesselEnhancementImageFilter: SetWStrength
    filter->SetIterations(atoi(argv[7]));
    filter->SetDarkObjectLightBackground(atoi(argv[8]));
    // filter->SetRecalculateVesselness();
    // filter->SetSensitivity();
    // filter->SetTimeStep();

    //// AnisotropicDiffusionVesselEnhancementImageFilter: SetSigmaMin, SetSigmaMax, SetNumberOfSigmaSteps
    std::vector<typename FilterType::Precision> scales;
    for (int i=0; i < NumberOfSigmaSteps; i++)
	scales.push_back(SigmaMin + (SigmaMax - SigmaMin) / (NumberOfSigmaSteps - 1) * i);

    filter->SetScales(scales);
    filter->VerboseOn();

    FilterWatcher watcher1(filter);
    try{
        filter->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }


    const typename OutputImageType::Pointer& output= filter->GetOutput();

    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer);
    writer->SetFileName(argv[2]);
    writer->SetInput(output);
    writer->SetUseCompression(atoi(argv[3]));
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
int dispatch_D(size_t dimensionType, int argc, char *argv[]){
    int res= 0;
    switch (dimensionType){
    //// even though VesselEnhancingDiffusion3DImageFilter is based on the Hessian matrix (which can be computed for 2D) it directly depends on 3 eigen-values, so only 3D is supported!
    case 3:
        res= DoIt<InputComponentType, InputPixelType, 3>(argc, argv);
        break;
    default:
        std::cerr << "Error: Images of dimension " << dimensionType << " are not handled!" << std::endl;
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
/*
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
*/
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
    if ( argc != 9 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Output_Image"
                  << " compress"
                  << " sigmaMin"
                  << " sigmaMax"
                  << " sigmaSteps"
                  << " iterations"
                  << " darkObject"
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






