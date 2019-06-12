#include <itkImageFileReader.h>
#include <itkResampleImageFilter.h>
#include <itkGaussianInterpolateImageFunction.h>
#include <itkImageFileWriter.h>

using ImageType = itk::Image<float, 3>;
using ReaderType = itk::ImageFileReader<ImageType>;
using FilterType = itk::ResampleImageFilter<ImageType, ImageType>;
using TransformType = itk::IdentityTransform<double, 3>;
using InterpolatorType = itk::GaussianInterpolateImageFunction<ImageType>;
using WriterType = itk::ImageFileWriter<ImageType>;

int main(int argc, char** argv)
{
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(argv[1]);
    reader->UpdateOutputInformation();
    ImageType::Pointer input = reader->GetOutput();

    typename TransformType::Pointer transform = TransformType::New();
    transform->SetIdentity();
    typename ImageType::SpacingType outputSpacing;
    outputSpacing[0]= 0.5;
    outputSpacing[1]= 0.4;
    outputSpacing[2]= 0.3;

    const typename ImageType::SpacingType& inputSpacing= input->GetSpacing();
    const typename ImageType::SizeType& inputSize= input->GetLargestPossibleRegion().GetSize();
    typename ImageType::SizeType outputSize;

    typedef typename ImageType::SizeType::SizeValueType SizeValueType;
    for (unsigned int i= 0; i < 3; i++)
        outputSize[i]= static_cast<SizeValueType>((double) inputSize[i] * inputSpacing[i] / outputSpacing[i]);


    typename InterpolatorType::Pointer interpolator= InterpolatorType::New();
    std::cerr << "Using interpolator: " << interpolator->GetNameOfClass() << std::endl;
    typename InterpolatorType::ArrayType sigma;
    for (unsigned int i= 0; i < 3; i++)
	sigma[i]= 0.8; //as suggested in pub: http://www.insight-journal.org/browse/publication/705
    interpolator->SetSigma(sigma);
    interpolator->SetAlpha(3.0);
    std::cerr << "Sigma: " << interpolator->GetSigma() << " Alpha: " << interpolator->GetAlpha() << std::endl;

    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(input);
    filter->SetTransform(transform);
    filter->SetInterpolator(interpolator);
    filter->SetOutputSpacing(outputSpacing);
    filter->SetSize(outputSize);
    filter->SetOutputOrigin(input->GetOrigin());
    filter->SetOutputDirection(input->GetDirection());
    filter->SetDefaultPixelValue(0);
    filter->ReleaseDataFlagOn();

    int chunks= atoi(argv[2]);
    std::stringstream sss;
    sss.str(""); sss << "out" << chunks << ".mhd";
   
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput(filter->GetOutput());
    try
    {
        writer->SetFileName(sss.str().c_str());
        writer->SetNumberOfStreamDivisions(chunks);
        writer->Update();
    }
    catch (itk::ExceptionObject& exc)
    {
        std::cout << exc;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
