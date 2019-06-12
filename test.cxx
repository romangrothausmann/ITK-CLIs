#include <itkImageFileReader.h>
#include <itkResampleImageFilter.h>
#include <itkGaussianInterpolateImageFunction.h>
#include <itkImageFileWriter.h>

using ImageType = itk::Image<int, 3>;
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
    outputSpacing.Fill(1.0);

    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(reader->GetOutput());
    filter->SetInterpolator(InterpolatorType::New());
    filter->SetSize(input->GetLargestPossibleRegion().GetSize());
    filter->SetOutputOrigin(input->GetOrigin());

    filter->SetTransform(transform);
    filter->SetOutputSpacing(outputSpacing);
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
