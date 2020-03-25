////image filter to apply itkShiftScaleImageFilter to get desired mean and Std based on mean and Std of ROI
//based on http://www.itk.org/Wiki/ITK/Examples/Developer/Minipipeline

#ifndef __itkROIShiftScaleImageFilter_h
#define __itkROIShiftScaleImageFilter_h
 
#include <itkImageToImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkShiftScaleImageFilter.h>

namespace itk{
    /** \class ROIShiftScaleImageFilter
     * \brief applies itkShiftScaleImageFilter to get desired mean and Std based on mean and Std of ROI
     *
     * \ingroup ImageFilters
     */
    template<typename TInputImage, typename TOutputImage>
	class ROIShiftScaleImageFilter: public ImageToImageFilter<TInputImage, TOutputImage>{

	public:
	/** Standard class typedefs. */
	typedef ROIShiftScaleImageFilter             Self;
	typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
	typedef SmartPointer<Self>        Pointer;
 
	/** Method for creation through the object factory. */
	itkNewMacro(Self);
 
	/** Run-time type information (and related methods). */
	itkTypeMacro(ROIShiftScaleImageFilter, ImageToImageFilter);
 
	/** Typedef to describe the output and input image region types. */
	typedef typename TInputImage::RegionType  InputImageRegionType;

	/** ImageDimension enumeration */
	itkStaticConstMacro(InputImageDimension, unsigned int, TInputImage::ImageDimension);
	itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

	/** Set/Get the output image region.
	 *  If any of the ExtractionRegion.Size = 0 for any particular dimension dim,
	 *  we have to collapse dimension dim.  This means the output image will have
	 *  'c' dimensions less than the input image, where c = number of
	 *  ExtractionRegion.Size = 0. */
	//void SetROI(InputImageRegionType ROI);
	itkSetMacro(ROI, InputImageRegionType);
	itkGetMacro(ROI, InputImageRegionType);

	itkSetMacro(DesiredMean, double);
	itkGetMacro(DesiredMean, double);

	itkSetMacro(DesiredStd, double);
	itkGetMacro(DesiredStd, double);

	/**  Filters  */
	typedef ExtractImageFilter<TInputImage, TOutputImage> ExtractFilterType;
	typedef StatisticsImageFilter<TInputImage> StatImageFilterType;
	typedef ShiftScaleImageFilter<TInputImage, TOutputImage> SSImageFilterType;
  
	protected:
	ROIShiftScaleImageFilter();
	~ROIShiftScaleImageFilter(){}
	void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

	/** Does the real work. */
	virtual void GenerateData();

	InputImageRegionType m_ROI;
	double m_DesiredMean, m_DesiredStd;

	private:
	ROIShiftScaleImageFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented
 
	};
    } //namespace ITK


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkROIShiftScaleImageFilter.hxx"
#endif

#endif // __itkROIShiftScaleImageFilter_h

