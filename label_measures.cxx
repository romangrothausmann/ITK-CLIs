////program for various label image measures: extents, stats, V, S, lol
//01: based on template.cxx


#include <complex>
#include <iomanip>

#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkChangeInformationImageFilter.h>

#include <itkStatisticsImageFilter.h>

#include <itkLabelImageToShapeLabelMapFilter.h>
#include <itkNumericTraits.h>
#include <itkShapeLabelObject.h>
#include <itkLabelMap.h>

#include <itkLabelOverlapMeasuresImageFilter.h>



template<typename InputComponentType, typename InputPixelType, size_t Dimension>
int DoIt(int argc, char *argv[]){


    typedef itk::Image<InputPixelType, Dimension>  InputImageType;


    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(argv[1]);
    reader->UpdateOutputInformation();

    std::cout << "input region index: " << reader->GetOutput()->GetLargestPossibleRegion().GetIndex()
              << "  size: " <<  reader->GetOutput()->GetLargestPossibleRegion().GetSize()
              << std::endl;

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

    typename ReaderType::Pointer reader2 = ReaderType::New();

    reader2->SetFileName(argv[2]);
    reader2->ReleaseDataFlagOn();
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

    typedef itk::ChangeInformationImageFilter<InputImageType> FilterType0;
    typename FilterType0::Pointer filter0= FilterType0::New();
    filter0->SetInput(reader->GetOutput());
    filter0->SetReferenceImage(reader2->GetOutput());
    filter0->UseReferenceImageOn();
    filter0->ChangeAll();
    FilterWatcher watcher0(filter0);
    try{ 
        filter0->Update();
        }
    catch(itk::ExceptionObject &ex){ 
	std::cerr << ex << std::endl;
	return EXIT_FAILURE;
	}

    const typename InputImageType::Pointer& input= filter0->GetOutput();

    typedef itk::StatisticsImageFilter<InputImageType> FilterType;
    typename FilterType::Pointer stat= FilterType::New();
    stat->SetInput(input);
    FilterWatcher watcher1(stat);

    try{ 
        stat->Update();
        }
    catch(itk::ExceptionObject &ex){ 
	std::cerr << ex << std::endl;
	return EXIT_FAILURE;
	}


    std::cout << "Min: " << +stat->GetMinimum() << " Max: " << +stat->GetMaximum() << " Mean: " << +stat->GetMean() << " Std: " << +stat->GetSigma() << " Variance: " << +stat->GetVariance() << " Sum: " << +stat->GetSum() << std::endl; //+ promotes variable to a type printable as a number (e.g. for char)

    typedef itk::LabelImageToShapeLabelMapFilter<InputImageType> FilterType2;//default output is a LabelMap instanciated with InputImageType::PixelType
    typename FilterType2::Pointer filter2= FilterType2::New();
    filter2->SetInput(input);
    filter2->ReleaseDataFlagOn();
    bool cp= atoi(argv[3]);
    filter2->SetComputePerimeter(cp);
    filter2->SetBackgroundValue(itk::NumericTraits<typename FilterType2::OutputImagePixelType>::max());

    FilterWatcher watcher2(filter2);
    try{
        filter2->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }

    typedef typename FilterType2::OutputImageType LabelMapType;
    typedef typename FilterType2::OutputImageType::LabelObjectType LabelObjectType;
    typedef typename FilterType2::OutputImageType::LabelType LabelType;
    typename LabelMapType::Pointer labelMap = filter2->GetOutput();
    std::cout << std::setw(10)
	      << "Label" << "\t" << std::setw(10)
	      << "Voxel" << "\t" << std::setw(10)
	      << "Vphy";
    if (cp)
	std::cout << "\t" << std::setw(10) << "Aphy";
    
    std::cout << std::endl;

    const LabelObjectType* labelObject;
    for(LabelType label= 0; label < labelMap->GetNumberOfLabelObjects(); label++){//SizeValueType == LabelType //GetNthLabelObject starts with 0 and ends at GetNumberOfLabelObjects()-1!!!

        labelObject= labelMap->GetNthLabelObject(label);//GetNthLabelObject essential in case of not consecutive labels! (compare: http://www.itk.org/Doxygen47/html/classitk_1_1BinaryImageToShapeLabelMapFilter.html and http://www.itk.org/Doxygen47/html/classitk_1_1LabelMap.html)
        std::cout
	    << std::setprecision(2) << std::right << std::setw(10)
            << +labelObject->GetLabel() << "\t" << std::setw(10); //essential in case of not consecutive labels!
        std::cout
            //<< labelObject->GetFeretDiameter() << "\t" << std::setw(10)
            //<< labelObject->GetElongation() << "\t" << std::setw(10)
            //<< labelObject->GetRoundness() << "\t" << std::setw(10)
            //<< labelObject->GetFlatness() << "\t" << std::setw(10)
            //<< labelObject->GetEquivalentSphericalRadius() << "\t" << std::setw(10)
            //<< labelObject->GetEquivalentSphericalPerimeter() << "\t" << std::setw(10)
            << +labelObject->GetNumberOfPixels() << "\t" << std::setw(10)
            << +labelObject->GetPhysicalSize();
        if (cp)
            std::cout << "\t" << std::setw(10) << +labelObject->GetPerimeter(); //takes physical voxel size into account!!!
        std::cout << std::endl;
        }


    typedef itk::LabelOverlapMeasuresImageFilter<InputImageType> FilterType3;
    typename FilterType3::Pointer filter3= FilterType3::New();
    filter3->SetSourceImage(input);
    filter3->SetTargetImage(reader2->GetOutput());
    filter3->ReleaseDataFlagOn();

    FilterWatcher watcher3(filter3);
    try{
        filter3->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }


    std::cout << std::setw(10)
	      << "Label" << "\t" << std::setw(10)
	      << "Total" << "\t" << std::setw(10)
	      << "Union" << "\t" << std::setw(10) // Jaccard
	      << "Mean" << "\t" << std::setw(10) // Dice
	      << "Simi" << "\t" << std::setw(10)
	      << "False-" << "\t" << std::setw(10)
	      << "False+" << "\t" << std::setw(10)
	      << std::endl;

    std::cout << "all" << "\t" << std::setw(10) // label 0 for total (not bg)
	      << std::fixed << std::setprecision(5) << std::right
	      << filter3->GetTotalOverlap() << "\t" << std::setw(10)
	      << filter3->GetUnionOverlap() << "\t" << std::setw(10)
	      << filter3->GetMeanOverlap() << "\t" << std::setw(10)
	      << filter3->GetVolumeSimilarity() << "\t" << std::setw(10)
	      << filter3->GetFalseNegativeError() << "\t" << std::setw(10)
	      << filter3->GetFalsePositiveError() << "\t" << std::setw(10)
	      << std::endl;

    typename FilterType3::MapType labelMap3 = filter3->GetLabelSetMeasures();
    typename FilterType3::MapType::const_iterator it;
    for(it= labelMap3.begin(); it != labelMap3.end(); ++it){
	if((*it).first == 0)
	    continue;

	int label= (*it).first;

	std::cout << label << "\t" << std::setw(10)
		  << std::fixed << std::setprecision(5) << std::right
		  << filter3->GetTargetOverlap(label) << "\t" << std::setw(10)
		  << filter3->GetUnionOverlap(label) << "\t" << std::setw(10)
		  << filter3->GetMeanOverlap(label) << "\t" << std::setw(10)
		  << filter3->GetVolumeSimilarity(label) << "\t" << std::setw(10)
		  << filter3->GetFalseNegativeError(label) << "\t" << std::setw(10)
		  << filter3->GetFalsePositiveError(label) << "\t" << std::setw(10)
		  << std::endl;
	}

    return EXIT_SUCCESS;

    }


template<typename InputComponentType, typename InputPixelType>
int dispatch_D(size_t dimensionType, int argc, char *argv[]){
    int res= EXIT_FAILURE;
    switch (dimensionType){
    // case 1:
    //     res= DoIt<InputComponentType, InputPixelType, 1>(argc, argv);
    //     break;
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
    int res= EXIT_FAILURE;
    //http://www.itk.org/Doxygen45/html/classitk_1_1ImageIOBase.html#abd189f096c2a1b3ea559bc3e4849f658
    //http://www.itk.org/Doxygen45/html/itkImageIOBase_8h_source.html#l00099
    //IOPixelType:: UNKNOWNPIXELTYPE, SCALAR, RGB, RGBA, OFFSET, VECTOR, POINT, COVARIANTVECTOR, SYMMETRICSECONDRANKTENSOR, DIFFUSIONTENSOR3D, COMPLEX, FIXEDARRAY, MATRIX

    switch (pixelType){
    case itk::ImageIOBase::SCALAR:{
        typedef InputComponentType InputPixelType;
        res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
        } break;
    // case itk::ImageIOBase::RGB:{
    //     typedef itk::RGBPixel<InputComponentType> InputPixelType;
    //     res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::RGBA:{
    //     typedef itk::RGBAPixel<InputComponentType> InputPixelType;
    //     res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::COMPLEX:{
    //     typedef std::complex<InputComponentType> InputPixelType;
    //     res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::VECTOR:{
    //     typedef itk::VariableLengthVector<InputComponentType> InputPixelType;
    //     res= dispatch_D<InputComponentType, InputPixelType>(dimensionType, argc, argv);
    //     } break;
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
    // case itk::ImageIOBase::FLOAT:{        // float32
    //     typedef float InputComponentType;
    //     res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
    //     } break;
    // case itk::ImageIOBase::DOUBLE:{       // float64
    //     typedef double InputComponentType;
    //     res= dispatch_pT<InputComponentType>(pixelType, dimensionType, argc, argv);
    //     } break;
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
    if ( argc != 4 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Comp_Image"
		  << " <bool: calculate perimeter>"
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






