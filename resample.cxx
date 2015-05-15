////program for itkResampleImageFilter
//01: based on template.cxx and http://itk.org/Wiki/ITK/Examples/ImageProcessing/ResampleSegmentedImage


#include <complex>

#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
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

    if( argc != 4 + 1*Dimension + 1){
	fprintf(stderr, "4 + 1*Dimension = %d parameters are needed!\n", 4 + 1*Dimension);
	return EXIT_FAILURE;
	}
	
    typedef InputPixelType  OutputPixelType;
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

    typename InputImageType::Pointer input= reader->GetOutput();


    typedef itk::IdentityTransform<TCoordRep, Dimension> TransformType;
    typename TransformType::Pointer transform= TransformType::New();
    transform->SetIdentity();

    const typename InputImageType::SpacingType& inputSpacing= reader->GetOutput()->GetSpacing();

    typename InputImageType::SpacingType outputSpacing;
    for (unsigned int i= 0; i < Dimension; i++)
        outputSpacing[i]= atof(argv[5+i]);

    const typename InputImageType::SizeType& inputSize= reader->GetOutput()->GetLargestPossibleRegion().GetSize();
    typename OutputImageType::SizeType outputSize;

    for (unsigned int i= 0; i < Dimension; i++)
	outputSize[i]= (double) inputSize[i] * inputSpacing[i] / outputSpacing[i];


    typedef itk::ResampleImageFilter<InputImageType, OutputImageType> FilterType;
    typename FilterType::Pointer filter= FilterType::New();
    filter->SetInput(input);
    filter->SetTransform(transform);
    filter->SetInterpolator(interpolator);
    filter->SetOutputSpacing(outputSpacing);
    filter->SetSize(outputSize);
    filter->SetDefaultPixelValue(itk::NumericTraits<InputPixelType>::Zero);
    filter->ReleaseDataFlagOn();
    //filter->InPlaceOn();//not available

    FilterWatcher watcher1(filter);
    try{
        filter->Update();
        }
    catch(itk::ExceptionObject &ex){
	std::cerr << ex << std::endl;
	return EXIT_FAILURE;
	}


    typename OutputImageType::Pointer output= filter->GetOutput();

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

template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){
    int res= 0;

    typedef double TCoordRep;
    typedef double TCoefficientType;
    typedef itk::Image<InputPixelType, Dimension>  InputImageType;

    switch(atoi(argv[4])){
    case 0:{
	typedef itk::NearestNeighborInterpolateImageFunction<InputImageType, TCoordRep> InterpolatorType;
    	typename InterpolatorType::Pointer interpolator= InterpolatorType::New();
    	std::cerr << "Using interpolator: " << interpolator->GetNameOfClass() << std::endl;
    	res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, TCoordRep, InterpolatorType>(argc, argv, interpolator);
    	}break;
    case 1:{
	typedef itk::LinearInterpolateImageFunction<InputImageType, TCoordRep> InterpolatorType;
    	typename InterpolatorType::Pointer interpolator= InterpolatorType::New();
    	std::cerr << "Using interpolator: " << interpolator->GetNameOfClass() << std::endl;
    	res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, TCoordRep, InterpolatorType>(argc, argv, interpolator);
    	}break;
    case 2:{
	typedef itk::BSplineInterpolateImageFunction<InputImageType, TCoordRep, TCoefficientType> InterpolatorType;
    	typename InterpolatorType::Pointer interpolator= InterpolatorType::New();
    	std::cerr << "Using interpolator: " << interpolator->GetNameOfClass() << std::endl;
    	interpolator->SetSplineOrder(3);
    	std::cerr << "Spline order: " << interpolator->GetSplineOrder() << std::endl;
    	res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, TCoordRep, InterpolatorType>(argc, argv, interpolator);
    	}break;
    case 3:{
	typedef itk::GaussianInterpolateImageFunction<InputImageType, TCoordRep> InterpolatorType;
    	typename InterpolatorType::Pointer interpolator= InterpolatorType::New();
    	std::cerr << "Using interpolator: " << interpolator->GetNameOfClass() << std::endl;
	interpolator->SetSigma(1.0);
	interpolator->SetAlpha(3.0);
	fprintf(stderr, "Sigma: %f;  Alpha: %f\n", interpolator->GetSigma(), interpolator->GetAlpha());
    	res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, TCoordRep, InterpolatorType>(argc, argv, interpolator);
    	}break;
    case 4:{
	typedef itk::LabelImageGaussianInterpolateImageFunction<InputImageType, TCoordRep> InterpolatorType;
    	typename InterpolatorType::Pointer interpolator= InterpolatorType::New();
    	std::cerr << "Using interpolator: " << interpolator->GetNameOfClass() << std::endl;
	interpolator->SetSigma(1.0);//as in: http://itk.org/Wiki/ITK/Examples/ImageProcessing/ResampleSegmentedImage
	interpolator->SetAlpha(3.0);
	fprintf(stderr, "Sigma: %f;  Alpha: %f\n", interpolator->GetSigma(), interpolator->GetAlpha());
    	res= DoIt2<InputComponentType, InputPixelType, Dimension, InputImageType, TCoordRep, InterpolatorType>(argc, argv, interpolator);
    	}break;
    case 5:{//from: http://www.itk.org/Doxygen/html/Examples_2Filtering_2ResampleImageFilter8_8cxx-example.html#_a7
	typedef itk::ConstantBoundaryCondition<InputImageType> BoundaryConditionType;
	const unsigned int WindowRadius = 5;
	typedef itk::Function::HammingWindowFunction<WindowRadius> WindowFunctionType;
    	typedef itk::WindowedSincInterpolateImageFunction<InputImageType, WindowRadius, WindowFunctionType, BoundaryConditionType, TCoordRep> InterpolatorType;
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
  // case itk::ImageIOBase::RGB:{ //does not work with: BSplineInterpolateImageFunction
  //   typedef itk::RGBPixel<InputComponentType> InputPixelType;
  //   res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
  // } break;
  // case itk::ImageIOBase::RGBA:{
  //   typedef itk::RGBAPixel<InputComponentType> InputPixelType;
  //   res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
  // } break;
  // case itk::ImageIOBase::COMPLEX:{ //does not work with: NearestNeighborInterpolateImageFunction
  //   typedef std::complex<InputComponentType> InputPixelType;
  //   res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
  // } break;
  // case itk::ImageIOBase::VECTOR:{ //does not work with: NearestNeighborInterpolateImageFunction
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
		  << " compress"
		  << " Interpolator_Type"
		  << " spacing..."
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






