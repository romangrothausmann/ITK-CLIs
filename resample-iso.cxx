////program to resample to iso voxels with additional smoothing of subsampled region
//01: based on resample.cxx and Examples/Filtering/ResampleVolumesToBeIsotropic.cxx
//02: single program for SDI and noSDI (based on resample.cxx)


#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkCastImageFilter.h>
#include <itkRecursiveGaussianImageFilter.h> // can stream: https://github.com/InsightSoftwareConsortium/ITK/blob/d5158e04f00cc61249cc1043147bc9a99f055d99/Modules/Filtering/ImageFilterBase/include/itkRecursiveSeparableImageFilter.hxx#L189-L193
#include <itkIdentityTransform.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkGaussianInterpolateImageFunction.h>
#include <itkLabelImageGaussianInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkWindowedSincInterpolateImageFunction.h>
#include <itkConstantBoundaryCondition.h>
#include <itkResampleImageFilter.h>
#include <itkImageFileWriter.h>




template<typename InputComponentType, typename InputPixelType, size_t Dimension, typename InputImageType, typename TCoordRep, typename InterpolatorType>
int DoIt2(int argc, char *argv[], InterpolatorType* interpolator){

    typedef float           InternalPixelType;
    typedef InputPixelType  OutputPixelType;
    typedef itk::Image<InternalPixelType, Dimension> InternalImageType;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;


    int CompChunk= atoi(argv[3]);
    bool noSDI= CompChunk <= 1; // SDI only if CompChunk > 1

    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(argv[1]);
    reader->ReleaseDataFlagOn();
    if(noSDI){
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
	}
    else{
	reader->UpdateOutputInformation();
	}

    const typename InputImageType::Pointer& input= reader->GetOutput();

    ////store header data before reader releases its output after it got used by the cast filter
    const typename InputImageType::SpacingType& inputSpacing= input->GetSpacing();
    const typename InputImageType::PointType& inputOrigin= input->GetOrigin();
    const typename InputImageType::DirectionType& inputDirection= input->GetDirection();
    const typename InputImageType::SizeType& inputSize= input->GetLargestPossibleRegion().GetSize();

    double isoSpacing;
    if(argc > 5)
	isoSpacing= atof(argv[5]);
    else{
	isoSpacing = std::sqrt(inputSpacing[0] * inputSpacing[Dimension-1]);
	std::cerr << "Using geometric mean between first and last spacing: " << isoSpacing << std::endl;
	}
    for (unsigned int i= 0; i < Dimension - 1; i++)
	if(isoSpacing < inputSpacing[i]){
	    fprintf(stderr, "Filter only intended for down-smapling but specified spacing (%e) is smaller than input spacing (%e) in dimension %d.\n", isoSpacing, inputSpacing[i], i);
	    return EXIT_FAILURE;
	    }

    double sigmaX= isoSpacing - inputSpacing[0]; // using diff, see e.g. https://discourse.itk.org/t/resampling-to-isotropic-signal-processing-theory/1403/13
    double sigmaY= isoSpacing - inputSpacing[1]; // using diff, see e.g. https://discourse.itk.org/t/resampling-to-isotropic-signal-processing-theory/1403/13

    std::cerr << "Using sigmaX (physical units): " << sigmaX << std::endl;
    std::cerr << "Using sigmaY (physical units): " << sigmaY << std::endl;

    typename InputImageType::SpacingType outputSpacing;
    for (unsigned int i= 0; i < Dimension; i++)
        outputSpacing[i]= isoSpacing;

    typename OutputImageType::SizeType outputSize;
    typedef typename InputImageType::SizeType::SizeValueType SizeValueType;
    for (unsigned int i= 0; i < Dimension; i++)
        outputSize[i]= static_cast<SizeValueType>((double) inputSize[i] * inputSpacing[i] / outputSpacing[i]);


    typedef itk::CastImageFilter<InputImageType, InternalImageType> CastFilterType;
    typename CastFilterType::Pointer  caster=  CastFilterType::New();
    caster->SetInput(input);
    caster->ReleaseDataFlagOn();
    caster->InPlaceOn();
    if(noSDI)
	FilterWatcher watcherC(caster);

    typedef itk::RecursiveGaussianImageFilter<InternalImageType, InternalImageType> GaussianFilterType;

    typename GaussianFilterType::Pointer smootherX = GaussianFilterType::New();
    smootherX->SetInput(caster->GetOutput());
    smootherX->SetSigma(sigmaX); // in physical units
    smootherX->SetDirection(0);
    smootherX->ReleaseDataFlagOn();
    smootherX->InPlaceOn();
    if(noSDI)
	FilterWatcher watcherX(smootherX);

    typename GaussianFilterType::Pointer smootherY = GaussianFilterType::New();
    smootherY->SetInput(smootherX->GetOutput());
    smootherY->SetSigma(sigmaY); // in physical units
    smootherY->SetDirection(1);
    smootherY->ReleaseDataFlagOn();
    smootherY->InPlaceOn();
    if(noSDI)
	FilterWatcher watcherY(smootherY);


    typedef itk::IdentityTransform<TCoordRep, Dimension> TransformType;
    typename TransformType::Pointer transform= TransformType::New();
    transform->SetIdentity();

    typedef itk::ResampleImageFilter<InternalImageType, OutputImageType> FilterType;
    typename FilterType::Pointer filter= FilterType::New();
    if(atoi(argv[4]))//expecting 0 <=> NN
	filter->SetInput(smootherY->GetOutput());
    else
	filter->SetInput(caster->GetOutput());
    filter->SetTransform(transform);
    filter->SetInterpolator(interpolator);
    filter->SetOutputSpacing(outputSpacing);
    filter->SetSize(outputSize);
    filter->SetOutputOrigin(inputOrigin);//essential for images created with e.g. itkExtractImageFilter
    filter->SetOutputDirection(inputDirection);
    filter->SetDefaultPixelValue(itk::NumericTraits<InputPixelType>::Zero);
    filter->ReleaseDataFlagOn();

    if(noSDI){
	FilterWatcher watcher1(filter);
	try{
	    filter->Update();
	    }
	catch(itk::ExceptionObject &ex){
	    std::cerr << ex << std::endl;
	    return EXIT_FAILURE;
	    }
	}


    const typename OutputImageType::Pointer& output= filter->GetOutput();

    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer);
    writer->SetFileName(argv[2]);
    writer->SetInput(output);
    if(noSDI){
	writer->SetUseCompression(CompChunk);
	}
    else{
	writer->UseCompressionOff(); // writing compressed is not supported when streaming!
	writer->SetNumberOfStreamDivisions(CompChunk);
	}
    try{
        writer->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    return EXIT_SUCCESS;

    }


template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){
    int res= 0;

    typedef double TCoordRep;
    typedef double TCoefficientType;
    typedef float           InternalPixelType;
    typedef itk::Image<InternalPixelType, Dimension> InternalImageType;
    typedef itk::Image<InputPixelType, Dimension>  InputImageType;

    int opt= atoi(argv[4]);
    switch(opt){
    case 0:{
        typedef itk::NearestNeighborInterpolateImageFunction<InternalImageType, TCoordRep> InterpolatorType;
        typename InterpolatorType::Pointer interpolator= InterpolatorType::New();
        std::cerr << "Using interpolator: " << interpolator->GetNameOfClass() << std::endl;
        res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, TCoordRep, InterpolatorType>(argc, argv, interpolator);
        }break;
    case 1:{
        typedef itk::LinearInterpolateImageFunction<InternalImageType, TCoordRep> InterpolatorType;
        typename InterpolatorType::Pointer interpolator= InterpolatorType::New();
        std::cerr << "Using interpolator: " << interpolator->GetNameOfClass() << std::endl;
        res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, TCoordRep, InterpolatorType>(argc, argv, interpolator);
        }break;
    case 2 ... 5: { // https://stackoverflow.com/questions/4494170/grouping-switch-statement-cases-together#28292802
        typedef itk::BSplineInterpolateImageFunction<InternalImageType, TCoordRep, TCoefficientType> InterpolatorType;
        typename InterpolatorType::Pointer interpolator= InterpolatorType::New();
        std::cerr << "Using interpolator: " << interpolator->GetNameOfClass() << std::endl;
        interpolator->SetSplineOrder(opt);
        std::cerr << "Spline order: " << interpolator->GetSplineOrder() << std::endl;
        res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, TCoordRep, InterpolatorType>(argc, argv, interpolator);
        }break;
    case 10:{
        typedef itk::GaussianInterpolateImageFunction<InternalImageType, TCoordRep> InterpolatorType;
        typename InterpolatorType::Pointer interpolator= InterpolatorType::New();
        std::cerr << "Using interpolator: " << interpolator->GetNameOfClass() << std::endl;
        typename InterpolatorType::ArrayType sigma;
        for (unsigned int i= 0; i < Dimension; i++)
            sigma[i]= 0.8; //as suggested in pub: http://www.insight-journal.org/browse/publication/705
        interpolator->SetSigma(sigma);
        interpolator->SetAlpha(3.0);
        std::cerr << "Sigma: " << interpolator->GetSigma() << " Alpha: " << interpolator->GetAlpha() << std::endl;
        res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, TCoordRep, InterpolatorType>(argc, argv, interpolator);
        }break;
    case 11:{
        typedef itk::LabelImageGaussianInterpolateImageFunction<InternalImageType, TCoordRep> InterpolatorType;
        typename InterpolatorType::Pointer interpolator= InterpolatorType::New();
        std::cerr << "Using interpolator: " << interpolator->GetNameOfClass() << std::endl;
        typename InterpolatorType::ArrayType sigma;
        for (unsigned int i= 0; i < Dimension; i++)
            sigma[i]= 0.8; //as suggested in pub: http://www.insight-journal.org/browse/publication/705
        interpolator->SetSigma(sigma);
        interpolator->SetAlpha(3.0);
        std::cerr << "Sigma: " << interpolator->GetSigma() << " Alpha: " << interpolator->GetAlpha() << std::endl;
        res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, TCoordRep, InterpolatorType>(argc, argv, interpolator);
        }break;
    case 20:{//from: http://www.itk.org/Doxygen/html/Examples_2Filtering_2ResampleImageFilter8_8cxx-example.html#_a7
        typedef itk::ConstantBoundaryCondition<InternalImageType> BoundaryConditionType;
        const unsigned int WindowRadius = 5;
        typedef itk::Function::HammingWindowFunction<WindowRadius> WindowFunctionType;
        typedef itk::WindowedSincInterpolateImageFunction<InternalImageType, WindowRadius, WindowFunctionType, BoundaryConditionType, TCoordRep> InterpolatorType;
        typename InterpolatorType::Pointer interpolator= InterpolatorType::New();
        std::cerr << "Using interpolator: " << interpolator->GetNameOfClass() << std::endl;
        fprintf(stderr, "With a %s and a window size of: %d\n", "HammingWindowFunction", WindowRadius);//no GetNameOfClass() for itkWindowFunction: http://public.kitware.com/pipermail/insight-users/2004-July/009440.html
        res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, TCoordRep, InterpolatorType>(argc, argv, interpolator);
        }break;
    default:
        std::cerr << "unknown interpolation type." << std::endl;
        res= EXIT_FAILURE;
        break;
        }//switch
    return res;
    }

template<typename InputComponentType, typename InputPixelType>
int dispatch_D(size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    switch (dimensionType){
    // case 2:
    //     res= DoIt<InputComponentType, InputPixelType, 2>(argc, argv);
    //     break;
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
    int res= EXIT_FAILURE;
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
    int res= EXIT_FAILURE;

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
    if ( argc < 5 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Output_Image"
                  << " compress|stream-chunks"
                  << " Interpolator_Type"
                  << " [iso-spacing]"
                  << std::endl;

        std::cerr << std::endl;
        std::cerr << " no-compress: 0, compress: 1, stream > 1" << std::endl;
        return EXIT_FAILURE;
        }

    int CompChunk= atoi(argv[3]);
    std::cerr << std::endl;
    if(CompChunk == 0){
	std::cerr << "Employing no compression and no streaming." << std::endl;
	}
    else if (CompChunk == 1){
	std::cerr << "Employing compression (streaming not possible then)." << std::endl;
	}
    else if (CompChunk > 1){
	std::cerr << "Employing streaming (compression not possible then)." << std::endl;
	}
    else {
	std::cerr << "compress|stream-chunks must be a positive integer" << std::endl;
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






