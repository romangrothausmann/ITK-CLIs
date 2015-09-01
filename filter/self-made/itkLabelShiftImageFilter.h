////image filter to shift labels to last the last maximum of a previous call
//based on http://www.itk.org/Wiki/ITK/Examples/Developer/Minipipeline

#ifndef __itkLabelShiftImageFilter_h
#define __itkLabelShiftImageFilter_h
 
#include <itkImageToImageFilter.h>
#include <itkAddImageFilter.h>

namespace itk{
    /** \class LabelShiftImageFilter
     * \brief image filter to shift labels to last the last maximum of a previous call
     *
     * \ingroup ImageFilters
     */
    template<typename TInputImage, typename TOutputImage>
	class LabelShiftImageFilter: public ImageToImageFilter<TInputImage, TOutputImage>{

	public:
	/** Standard class typedefs. */
	typedef LabelShiftImageFilter             Self;
	typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
	typedef SmartPointer<Self>        Pointer;
 
	/** Method for creation through the object factory. */
	itkNewMacro(Self);
 
	/** Run-time type information (and related methods). */
	itkTypeMacro(LabelShiftImageFilter, ImageToImageFilter);
 
	/** Typedef to describe the output and input image region types. */
	typedef typename TInputImage::RegionType  InputImageRegionType;
	typedef typename TOutputImage::RegionType OutputImageRegionType;

	/** ImageDimension enumeration */
	itkStaticConstMacro(InputImageDimension, unsigned int, TInputImage::ImageDimension);
	itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

	itkGetMacro(LastMax, typename TInputImage::PixelType);

	/**  Filters  */
	typedef AddImageFilter<TInputImage, TOutputImage> AddFilterType;
  
	protected:
	LabelShiftImageFilter();
	~LabelShiftImageFilter(){}
	void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

	/** Does the real work. */
	virtual void GenerateData();
	//virtual void ThreadedGenerateData(const OutputImageRegionType &, ThreadIdType);

	typename TInputImage::PixelType m_LastMax;

	typename AddFilterType::Pointer m_AddImageFilter;

	private:
	LabelShiftImageFilter(const Self &); //purposely not implemented
	void operator=(const Self &);  //purposely not implemented
 
	};
    } //namespace ITK


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLabelShiftImageFilter.hxx"
#endif

#endif // __itkLabelShiftImageFilter_h

