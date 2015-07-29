////program to segment a binray image along a minimal path
//01: based on min-path_fm.cxx and skeletonize_3D.cxx

////ToDo:


#include <string>
#include <sstream>

#include "itkFilterWatcher.h"
#include <itkImageFileReader.h>
#include <itkParabolicOpenImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkSpeedFunctionToPathFilter.h>
#include <itkMorphologicalDistanceTransformImageFilter.h>
#include <itkIterateNeighborhoodOptimizer.h>
#include <itkGradientDescentOptimizer.h>
#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkPathIterator.h>
#include <itkPolyLineParametricPath.h>
#include <itkParabolicDilateImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkImageFileWriter.h>

#include <itkMesh.h>
#include <itkLineCell.h>
#include <itkMeshFileWriter.h>
#include <itkVTKPolyDataMeshIO.h>
#include <itkMetaDataObject.h>


////this class just creates an alias from GetCurrentValue() to GetValue() in order to be API compliant
namespace itk{
    class MyIterateNeighborhoodOptimizer: public IterateNeighborhoodOptimizer{
    public:

	/** Standard class typedefs. */
	typedef MyIterateNeighborhoodOptimizer   Self;
	typedef IterateNeighborhoodOptimizer Superclass;
	typedef SmartPointer<Self>             Pointer;
	typedef SmartPointer<const Self>       ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro(Self); //essential for typedef creation, needs all typedefs above!

	/** Run-time type information (and related methods). */
	itkTypeMacro(MyIterateNeighborhoodOptimizer, IterateNeighborhoodOptimizer);

	// "rename" GetCurrentValue() to GetValue()
	double GetValue(){this->GetCurrentValue();};
	};
    }


template<typename OptimizerType>
void FilterEventHandlerITK(itk::Object *caller, const itk::EventObject &event, void*){

    const itk::ProcessObject* filter = static_cast<const itk::ProcessObject*>(caller);

    if(itk::ProgressEvent().CheckEvent(&event))
        fprintf(stderr, "\r%s progress: %5.1f%%", filter->GetNameOfClass(), 100.0 * filter->GetProgress());//stderr is flushed directly
    else if(itk::IterationEvent().CheckEvent(&event)){
        OptimizerType* optimizer= dynamic_cast<OptimizerType *>(caller);
	fprintf(stderr, "\r%5d: ", optimizer->GetCurrentIteration());
	if(optimizer->GetValue() < 1e18){
	    fprintf(stderr, "%7.3f ", optimizer->GetValue());
	    std::cerr << optimizer->GetCurrentPosition();
	    }
        }
    else if(itk::EndEvent().CheckEvent(&event))
        std::cerr << std::endl;
    }


template<typename InputComponentType, typename InputPixelType, size_t Dimension, typename OptimizerType>
int DoIt2(int argc, char *argv[], OptimizerType* optimizer){

    const char offset= 8;
    if((argc - offset) % Dimension){
        fprintf(stderr, "%d + n*Dimension  parameters are needed!\n", offset-1);
        return EXIT_FAILURE;
        }

    char* outPrefix= argv[2];
    std::stringstream sss;

#ifdef USE_FLOAT
    typedef float   SpeedPixelType;
    typedef float   DMPixelType;
    std::cerr << "Using single precision (float)." << std::endl;
#else
    typedef double   SpeedPixelType;
    typedef double   DMPixelType;
    std::cerr << "Using double precision (double)." << std::endl;
#endif

    typedef uint8_t  OutputPixelType;

    typedef itk::Image<InputPixelType, Dimension>  InputImageType;
    typedef itk::Image<SpeedPixelType, Dimension>  SpeedImageType;
    typedef itk::Image<DMPixelType, Dimension>  DMImageType;
    typedef itk::Image<OutputPixelType, Dimension>  OutputImageType;


    typedef itk::ImageFileReader<InputImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName(argv[1]);
    reader->ReleaseDataFlagOn();//save mem, reader will re-execute later on
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

    ////check for isotropic image spacing as dm later on uses this image spacing
    for(int i= 0; i < Dimension; i++)
	if(reader->GetOutput()->GetSpacing()[0] != reader->GetOutput()->GetSpacing()[i]){
	    std::cerr << "Image spacing not isotropic! This has not been tested, aborting." << std::endl;
	    return EXIT_FAILURE;
	    }
		

    typename SpeedImageType::Pointer speed;
	{////scoped to save mem and for consitency
	// scale image values to be in [0; 1]
	typedef typename itk::RescaleIntensityImageFilter<InputImageType, SpeedImageType> RescaleFilterType;
	typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput(reader->GetOutput());
	rescaleFilter->SetOutputMinimum(0.0);
	rescaleFilter->SetOutputMaximum(1.0);
	rescaleFilter->ReleaseDataFlagOn();//save mem
	FilterWatcher watcherR(rescaleFilter);

	if(atof(argv[4]) > 0){
	    std::cerr << "Sigma > 0, expecting binary image as input, i.e. using ParabolicOpenImageFilter to create speed map." << std::endl;

	    typedef itk::ParabolicOpenImageFilter<SpeedImageType, SpeedImageType> SmoothFilterType;
	    typename SmoothFilterType::Pointer smoother= SmoothFilterType::New();
	    smoother->SetInput(rescaleFilter->GetOutput());
	    smoother->SetScale(atof(argv[4]));
	    smoother->UseImageSpacingOff();
	    smoother->ReleaseDataFlagOn();//save mem
	    FilterWatcher watcherSM(smoother);

	    try {
		smoother->Update();
		}
	    catch (itk::ExceptionObject &ex){
		std::cerr << ex << std::endl;
		return EXIT_FAILURE;
		}

	    speed= smoother->GetOutput();
	    speed->DisconnectPipeline();
	    }
	else{
	    std::cerr << "Sigma <= 0, expecting speed map as input." << std::endl;
	    try {
		rescaleFilter->Update();
		}
	    catch (itk::ExceptionObject &ex){
		std::cerr << ex << std::endl;
		return EXIT_FAILURE;
		}

	    speed= rescaleFilter->GetOutput();
	    speed->DisconnectPipeline();
	    }
	}

    typedef itk::PolyLineParametricPath<Dimension> PathType;
    typedef itk::SpeedFunctionToPathFilter<SpeedImageType, PathType> PathFilterType;
    typedef typename PathFilterType::CostFunctionType::CoordRepType CoordRepType;

    // Create interpolator
    typedef itk::LinearInterpolateImageFunction<SpeedImageType, CoordRepType> InterpolatorType;
    typename InterpolatorType::Pointer interp = InterpolatorType::New();

    // Create cost function
    typename PathFilterType::CostFunctionType::Pointer cost = PathFilterType::CostFunctionType::New();
    cost->SetInterpolator(interp);

    itk::CStyleCommand::Pointer eventCallbackITK;
    eventCallbackITK = itk::CStyleCommand::New();
    eventCallbackITK->SetCallback(FilterEventHandlerITK<OptimizerType>);
    optimizer->AddObserver(itk::AnyEvent(), eventCallbackITK);

    // Create path filter
    typename PathFilterType::Pointer pathFilter = PathFilterType::New();
    pathFilter->SetInput(speed); //needs the image values to be scaled to [0; 1]
    pathFilter->SetCostFunction(cost);
    pathFilter->SetOptimizer(optimizer);
    pathFilter->SetTerminationValue(2.0);

    // Setup path points
    typename SpeedImageType::IndexType start, end;
    typename PathFilterType::PointType startP, endP;
    typename PathFilterType::PathInfo info;

    for(int j= 0; j < Dimension; j++){
        startP[j]= atof(argv[j+offset]+1);//+1 to skip prefix-letter
	start[j]= startP[j] - 1;//-1 for ITK-Snap compat
        }

    for(int j= 0; j < Dimension; j++){
        endP[j]= atof(argv[j+Dimension+offset]+1);//+1 to skip prefix-letter
 	end[j]= endP[j] - 1;//-1 for ITK-Snap compat
        }

    if(argv[offset][0]=='v')
	speed->TransformIndexToPhysicalPoint(start, startP);//overwrites startP
    if(argv[Dimension+offset][0]=='v')
	speed->TransformIndexToPhysicalPoint(end, endP);//overwrites endP

    info.SetStartPoint(startP);
    info.SetEndPoint(endP);
    std::cerr << "S: " << startP << std::endl;	
    std::cerr << "E: " << endP << std::endl;
	
    const typename SpeedImageType::RegionType region= speed->GetLargestPossibleRegion();
    if(!region.IsInside(start)){std::cerr << "Start point not inside image region. Aborting!" << std::endl; return EXIT_FAILURE;}
    if(!region.IsInside(end)){std::cerr << "End point not inside image region. Aborting!" << std::endl; return EXIT_FAILURE;}

    for(int i= offset + 2*Dimension; i < argc; i+= Dimension){
	typename SpeedImageType::IndexType way;
	typename PathFilterType::PointType wayP;

	for(int j= 0; j < Dimension; j++){
	    wayP[j]= atof(argv[i+j]+1);//+1 to skip prefix-letter
	    way[j]= wayP[j] - 1;//-1 for ITK-Snap compat
	    }

	if(argv[i][0]=='v')
	    speed->TransformIndexToPhysicalPoint(way, wayP);//overwrites wayP

	info.AddWayPoint(wayP);
	std::cerr << "W: " << wayP << std::endl;	
	if(!region.IsInside(way)){std::cerr << "Way point not inside image region. Aborting!" << std::endl; return EXIT_FAILURE;}
	}

    pathFilter->AddPathInfo(info);
    FilterWatcher watcher(pathFilter); //filter reports no progress so far
    try {
        pathFilter->Update();
        }
    catch (itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }


    //// create dm to get max-inscribed sphere radius for the point data in the vector/mesh output and for the max-sphere segementation
    typedef itk::MorphologicalDistanceTransformImageFilter<InputImageType, DMImageType> DMFilterType;
    typename DMFilterType::Pointer dm= DMFilterType::New();
    dm->SetInput(reader->GetOutput());
    //dm->ReleaseDataFlagOn();
    dm->SqrDistOn();
    dm->SetUseImageSpacing(true);//use image spacing to make "MaxInscrSphereRadius" fit to resampled volumes

    FilterWatcher watcherDM(dm);
    try{ 
        dm->Update();
        }
    catch(itk::ExceptionObject &ex){ 
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }
    const typename DMImageType::Pointer& dmap= dm->GetOutput();


    //// create mesh to save in a VTK-file
    typedef typename itk::Mesh<float, Dimension>  MeshType;
    typename MeshType::Pointer  mesh = MeshType::New();

    typename MeshType::PointType mP;
    for (unsigned int i=0; i < pathFilter->GetNumberOfOutputs(); i++){

        // Get the path, coords are stored as continous index
        typename PathType::Pointer path = pathFilter->GetOutput(i);
        const typename PathType::VertexListType *vertexList = path->GetVertexList();

        for(unsigned int k = 0; k < vertexList->Size(); k++){
            speed->TransformContinuousIndexToPhysicalPoint(vertexList->GetElement(k), mP);
            mesh->SetPoint(k, mP);

	    ////remeber to set scaling to voxel-size when visualizing point data with sphere glyphs in paraview
	    typename SpeedImageType::IndexType index;
	    index.CopyWithRound(vertexList->GetElement(k));
            mesh->SetPointData(k, std::sqrt(dmap->GetPixel(index)));
            }
        }
    // mesh->GetPointData()->SetObjectName("MaxInscrSphereRadius");
    // itk::MetaDataDictionary & metaDic= mesh->GetPointData()->GetMetaDataDictionary();
    // itk::EncapsulateMetaData<std::string>(metaDic, "pointScalarDataName", "MaxInscrSphereRadius");
    // mesh->GetPointData()->SetMetaDataDictionary(metaDic);

    std::cout << "# of mesh points: " << mesh->GetNumberOfPoints() << std::endl;

    //// create connecting lines
    typedef typename MeshType::CellType CellType;
    typedef typename itk::LineCell<CellType> LineType;
    typedef typename CellType::CellAutoPointer  CellAutoPointer;

    //// from: http://www.itk.org/Doxygen/html/Examples_2DataRepresentation_2Mesh_2Mesh3_8cxx-example.html
    if(mesh->GetNumberOfPoints() > 1){
        const unsigned int numberOfCells = mesh->GetNumberOfPoints() - 1;
        typename CellType::CellAutoPointer line;
        for(size_t cellId=0; cellId < numberOfCells; cellId++){
            line.TakeOwnership(new LineType);
            line->SetPointId(0, cellId);
            line->SetPointId(1, cellId+1);
            mesh->SetCell(cellId, line);
            }
        }

    std::cout << "# of mesh cells: " << mesh->GetNumberOfCells() << std::endl;

    typedef typename itk::MeshFileWriter<MeshType> MeshWriterType;
    typename MeshWriterType::Pointer mwriter = MeshWriterType::New();

    FilterWatcher watcherMO(mwriter);
    sss.str(""); sss << outPrefix << ".vtk"; //vtp not supported as of itk-4.8
    mwriter->SetFileName(sss.str().c_str());
    mwriter->SetInput(mesh);

    itk::VTKPolyDataMeshIO::Pointer mio= itk::VTKPolyDataMeshIO::New();
    itk::MetaDataDictionary & metaDic= mio->GetMetaDataDictionary();
    itk::EncapsulateMetaData<std::string>(metaDic, "pointScalarDataName", "MaxInscrSphereRadius");
    mwriter->SetMeshIO(mio);

    try{
	mwriter->Update();
        }
    catch(itk::ExceptionObject &ex){
        std::cerr << ex << std::endl;
        return EXIT_FAILURE;
        }


    // Allocate output image
    typename DMImageType::Pointer output = DMImageType::New();
    output->SetRegions(speed->GetLargestPossibleRegion());
    output->SetSpacing(speed->GetSpacing());
    output->SetOrigin(speed->GetOrigin());
    output->Allocate();
    output->FillBuffer(itk::NumericTraits<OutputPixelType>::Zero);


    // Rasterize path
    typedef itk::PathIterator<DMImageType, PathType> PathIteratorType;
    typedef itk::PathConstIterator<DMImageType, PathType> PathConstIteratorType;
    std::cout << "# of paths: " << pathFilter->GetNumberOfOutputs() << std::endl;
    for (unsigned int i= 0; i < pathFilter->GetNumberOfOutputs(); i++){

        // Get the path
        typename PathType::Pointer path = pathFilter->GetOutput(i);

        // Check path is valid
        if (path->GetVertexList()->Size() == 0){
            std::cerr << "WARNING: Path " << (i+1) << " contains no points!" << std::endl;
            continue;
            }

        printf("Path %3d contains %6d points.\n", i, path->GetVertexList()->Size());
        std::cerr << path->EndOfInput() << std::endl;

        // Iterate path and convert to image
        PathIteratorType it(output, path);
        PathConstIteratorType cit(dmap, path);
        int count= 0;
	typename PathConstIteratorType::PixelType value;
	typename PathConstIteratorType::PixelType min_v= itk::NumericTraits<typename PathConstIteratorType::PixelType>::max();
	typename PathConstIteratorType::PixelType max_v= itk::NumericTraits<typename PathConstIteratorType::PixelType>::min();
	typename PathConstIteratorType::IndexType min_i, max_i;
        for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++cit){//possibly more iterations than path->GetVertexList()->Size() !! can take long with large radii
	    value= cit.Get();
            it.Set(value);//mark center with squared hight, later used by ParabolicDilate
	    if(value < min_v){
		min_v= value;
		min_i= cit.GetIndex();
		}
	    if(value > max_v){
		max_v= value;
		max_i= cit.GetIndex();
		}
            }
	std::cerr << "Min radius along path: " << std::sqrt(min_v) << " [" << std::sqrt(min_v)/speed->GetSpacing()[0] << " v] @: " << min_i << std::endl;
	std::cerr << "Max radius along path: " << std::sqrt(max_v) << " [" << std::sqrt(max_v)/speed->GetSpacing()[0] << " v] @: " << max_i << std::endl;
	}
    dmap->ReleaseData();


    typedef itk::ParabolicDilateImageFilter<DMImageType, DMImageType> PDFilterType;
    typename PDFilterType::Pointer pd= PDFilterType::New();
    pd->SetInput(output);
    pd->SetScale(.5);//0.5: unscaled
    pd->SetUseImageSpacing(true);//needs to be concise with dm above
    pd->ReleaseDataFlagOn();
    FilterWatcher watcherPD(pd);

    typedef itk::BinaryThresholdImageFilter<DMImageType, OutputImageType> ThrFilterType;
    typename ThrFilterType::Pointer thr= ThrFilterType::New();
    thr->SetInput(pd->GetOutput());
    thr->SetUpperThreshold(0);
    thr->SetOutsideValue(1);
    thr->SetInsideValue(0);
    thr->ReleaseDataFlagOn();
    FilterWatcher watcherThr(thr);

    typedef itk::ImageFileWriter<OutputImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();

    FilterWatcher watcherO(writer);
    sss.str(""); sss << outPrefix << "_max-sphere.mha";
    writer->SetFileName(sss.str().c_str());
    writer->SetInput(thr->GetOutput());
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

    switch(atoi(argv[5])){
    case 0:{
	typedef itk::MyIterateNeighborhoodOptimizer OptimizerType;
        typename OptimizerType::Pointer optimizer= OptimizerType::New();
	optimizer->MinimizeOn();
	optimizer->FullyConnectedOn();

	typedef itk::Image<InputPixelType, Dimension>  InputImageType;
	typedef itk::ImageFileReader<InputImageType> ReaderType;
	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(argv[1]);
        reader->UpdateOutputInformation();

	typename OptimizerType::NeighborhoodSizeType size(Dimension);
	for (unsigned int i=0; i<Dimension; i++)
	    size[i] = reader->GetOutput()->GetSpacing()[i] * atof(argv[7]);
	optimizer->SetNeighborhoodSize(size);

        std::cerr << "Using interpolator: " << optimizer->GetNameOfClass() << " (ignoring iterations parameter)" << std::endl;
        res= DoIt2<InputComponentType, InputPixelType, Dimension, OptimizerType>(argc, argv, optimizer);
        }break;
    case 1:{
	typedef itk::GradientDescentOptimizer OptimizerType;
	typename OptimizerType::Pointer optimizer = OptimizerType::New();
	optimizer->SetNumberOfIterations(atoi(argv[6]));
	optimizer->SetLearningRate(atof(argv[7]));
        std::cerr << "Using interpolator: " << optimizer->GetNameOfClass() << std::endl;
        res= DoIt2<InputComponentType, InputPixelType, Dimension, OptimizerType>(argc, argv, optimizer);
        }break;
    case 2:{
	typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
	typename OptimizerType::Pointer optimizer = OptimizerType::New();
	optimizer->SetNumberOfIterations(atoi(argv[6]));
	optimizer->SetRelaxationFactor(.5);
	optimizer->SetMaximumStepLength(1.0);
	optimizer->SetMinimumStepLength(atof(argv[7]));
        std::cerr << "Using interpolator: " << optimizer->GetNameOfClass() << std::endl;
        res= DoIt2<InputComponentType, InputPixelType, Dimension, OptimizerType>(argc, argv, optimizer);
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
    int res= 0;
    switch (dimensionType){
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
    if ( argc < 8 ){
        std::cerr << "Missing Parameters: "
                  << argv[0]
                  << " Input_Image"
                  << " Output_Image_Base"
                  << " compress"
                  << " sigma"
                  << " optimizer"
                  << " iterations"
                  << " step-scale"
                  << " start-point..."
                  << " end-point..."
                  << " way-point..."
                  << std::endl;
        std::cerr << " Point coordinates are expected in voxel units (starting with 1)!" << std::endl;
        std::cerr << " step-scale will correspond to distance between points for a speed function ~1 along path" << std::endl;

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






