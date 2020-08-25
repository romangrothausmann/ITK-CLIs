////program to extract sub-images at cc-centers (cross-sections of a 3D skeleton, for stereological evaluations)
//01: based on template_2inputs.cxx and subimage_extract.cxx analyse_binary.cxx label_stack.cxx


#include <complex>

#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkExtractImageFilter.h>
#include <itkJoinSeriesImageFilter.h>
#include <itkBinaryImageToShapeLabelMapFilter.h>
#include <itkShapeLabelObject.h>
#include <itkLabelMap.h>
#include <itkImageFileWriter.h>



template<typename InputComponentType1, typename InputComponentType2, typename InputPixelType1, typename InputPixelType2, size_t Dimension>
int DoIt(int argc, char *argv[]){


    typedef InputPixelType1  OutputPixelType;

    typedef itk::Image<InputPixelType1, Dimension>  InputImageType1;
    typedef itk::Image<InputPixelType2, Dimension>  InputImageType2;
    typedef itk::Image<OutputPixelType, Dimension>  SliceImageType;
    typedef itk::Image<OutputPixelType, Dimension + 1>  OutputImageType;


    int CompChunk= 0;
    bool noSDI= true;

    typedef itk::ImageFileReader<InputImageType1> ReaderType1;
    typename ReaderType1::Pointer reader1 = ReaderType1::New();

    reader1->SetFileName(argv[1]);
    // do not release, will be needed multiple times: reader1->ReleaseDataFlagOn();
    if(noSDI){
	FilterWatcher watcherI1(reader1);
	watcherI1.QuietOn();
	watcherI1.ReportTimeOn();
	try{
	    reader1->Update();
	    }
	catch(itk::ExceptionObject &ex){
	    std::cerr << ex << std::endl;
	    return EXIT_FAILURE;
	    }
	}
    else{
	reader1->UpdateOutputInformation();
	}

    const typename InputImageType1::Pointer& input1= reader1->GetOutput();


    typedef itk::ImageFileReader<InputImageType2> ReaderType2;
    typename ReaderType2::Pointer reader2 = ReaderType2::New();

    reader2->SetFileName(argv[2]);
    reader2->ReleaseDataFlagOn();
    if(noSDI){
	FilterWatcher watcherI2(reader2);
	watcherI2.QuietOn();
	watcherI2.ReportTimeOn();
	try{
	    reader2->Update();
	    }
	catch(itk::ExceptionObject &ex){
	    std::cerr << ex << std::endl;
	    return EXIT_FAILURE;
	    }
	}
    else{
	reader2->UpdateOutputInformation();
	}

    const typename InputImageType2::Pointer& input2= reader2->GetOutput();

    unsigned int i;
    typename InputImageType1::IndexType desiredStart;
    typename InputImageType1::SizeType desiredSize;

    for (i= 0; i < Dimension; i++)
        desiredSize[i]= atoi(argv[4]);

    typedef itk::ExtractImageFilter<InputImageType1, SliceImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(input1);
    filter->SetDirectionCollapseToIdentity(); // This is required.
    FilterWatcher watcher0(filter);

    typedef itk::JoinSeriesImageFilter<SliceImageType, OutputImageType> JoinSeriesFilterType;
    typename JoinSeriesFilterType::Pointer joinSeries = JoinSeriesFilterType::New();
    joinSeries->SetOrigin(0); // origin for Dimension + 1
    joinSeries->SetSpacing(1); // spacing for Dimension + 1

    
    typedef itk::BinaryImageToShapeLabelMapFilter<InputImageType2> LMType;
    typename LMType::Pointer lm= LMType::New();
    lm->SetInput(reader2->GetOutput());
    // lm->SetInputForegroundValue();
    lm->FullyConnectedOn();
    lm->ComputePerimeterOff();
    lm->ReleaseDataFlagOn();

    if(noSDI){
	FilterWatcher watcher1(lm);
	try{
	    lm->Update();
	    }
	catch(itk::ExceptionObject &ex){
	    std::cerr << ex << std::endl;
	    return EXIT_FAILURE;
	    }
	}
    
    typedef typename LMType::OutputImageType LabelMapType;
    typedef typename LMType::OutputImageType::LabelObjectType LabelObjectType;
    typedef typename LMType::OutputImageType::LabelType LabelType;

    //// list of quantities see: http://www.itk.org/Doxygen/html/classitk_1_1ShapeLabelObject.html
    typename LabelMapType::Pointer labelMap = lm->GetOutput();

    const LabelObjectType* labelObject;
    for(LabelType label= 0; label < labelMap->GetNumberOfLabelObjects(); label++){//SizeValueType == LabelType //GetNthLabelObject starts with 0 and ends at GetNumberOfLabelObjects()-1!!!

        labelObject= labelMap->GetNthLabelObject(label);//using GetNthLabelObject to be save (even though the doc suggests otherwise (compare: http://www.itk.org/Doxygen47/html/classitk_1_1BinaryImageToShapeLabelMapFilter.html and http://www.itk.org/Doxygen47/html/classitk_1_1LabelMap.html)
	for (i= 0; i < Dimension; i++)
	    desiredStart[i]=  labelObject->GetCentroid()[i] - desiredSize[i] / 2;

	typename InputImageType1::RegionType desiredRegion(desiredStart, desiredSize);
	std::cerr << "desired region index: " << desiredRegion.GetIndex()
		  << "  size: " << desiredRegion.GetSize()
		  << std::endl;
	
	if(!reader1->GetOutput()->GetLargestPossibleRegion().IsInside(desiredRegion)){
	    std::cerr << "Desired region is not inside the largest possible input region, skipping connected component" << std::endl;
	    continue; // skip this region because it would need a boundary extension (e.g. padding, mirroring)
	    }

	filter->SetExtractionRegion(desiredRegion);

	try{
	    filter->Update();
	    }
	catch(itk::ExceptionObject &ex){
	    std::cerr << ex << std::endl;
	    continue;
	    }

	typename SliceImageType::Pointer subImg;
	subImg= filter->GetOutput();
	subImg->DisconnectPipeline();

	joinSeries->PushBackInput(subImg);
	}
	    
    
    const typename OutputImageType::Pointer& output= joinSeries->GetOutput();

    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer);
    writer->SetFileName(argv[3]);
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


template<typename InputComponentType1, typename InputComponentType2, typename InputPixelType1, typename InputPixelType2>
int dispatch_D(size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    switch (dimensionType){
    case 2:
        res= DoIt<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2, 2>(argc, argv);
        break;
    // case 3:
    //     res= DoIt<InputComponentType1, InputComponentType2, InputPixelType1, InputPixelType2, 3>(argc, argv);
    //     break;
    default:
        std::cerr << "Error: Images of dimension " << dimensionType << " are not handled!" << std::endl;
        break;
        }//switch
    return res;
    }

template<typename InputComponentType1, typename InputComponentType2, typename InputPixelType1>
int dispatch_pT2(itk::ImageIOBase::IOPixelType pixelType2, size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (pixelType2){
    case itk::ImageIOBase::SCALAR:{ // 1 component per pixel
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
    int res= EXIT_FAILURE;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (pixelType1){
    case itk::ImageIOBase::SCALAR:{ // 1 component per pixel
        typedef InputComponentType1 InputPixelType1;
        res= dispatch_pT2<InputComponentType1, InputComponentType2, InputPixelType1>(pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::COMPLEX:{ // 2 components per pixel
        typedef std::complex<InputComponentType1> InputPixelType1;
        res= dispatch_pT2<InputComponentType1, InputComponentType2, InputPixelType1>(pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::RGB:{ // 3 components per pixel, limited [0,1]
        typedef itk::RGBPixel<InputComponentType1> InputPixelType1;
        res= dispatch_pT2<InputComponentType1, InputComponentType2, InputPixelType1>(pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::RGBA:{ // 4 components per pixel, limited [0,1]
        typedef itk::RGBAPixel<InputComponentType1> InputPixelType1;
        res= dispatch_pT2<InputComponentType1, InputComponentType2, InputPixelType1>(pixelType2, dimensionType, argc, argv);
        } break;
    case itk::ImageIOBase::VECTOR:{
        typedef itk::VariableLengthVector<InputComponentType1> InputPixelType1;
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
    int res= EXIT_FAILURE;

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
    int res= EXIT_FAILURE;

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
    if ( argc != 5 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image1"
                  << " Input_Image2"
                  << " Output_Image"
                  << " crop-size"
                  << std::endl;

        std::cerr << std::endl;
        std::cerr << " no-compress: 0, compress: 1, stream > 1" << std::endl;
        return EXIT_FAILURE;
        }

    int CompChunk= atoi(argv[4]);
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






