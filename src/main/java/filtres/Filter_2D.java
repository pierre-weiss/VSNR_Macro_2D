package filtres;

import process.Double2DArray_2D;
import process.Double3DArray_2D;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;


public class Filter_2D {
	private FilterType_2D type;
	private Double3DArray_2D stackFilter;
	private Double2DArray_2D flt;
	private ImagePlus img;
	private ImageStack stack;
	private int p;
	private double alpha;
	private int width;
	private int height;
	private int slices;

	public Filter_2D(ImagePlus imgF, FilterType_2D tF, int widthF, int heightF) {
		this.img = IJ.createImage("Filterzzz", "32-bit", width, height, 1);
		this.img = new ImagePlus("Filter", imgF.getProcessor());
		this.width = widthF;
		this.height = heightF;
		this.type = tF;
		this.p = 2;
		this.alpha = 1;
		if (tF == FilterType_2D.DIRAC) {
			this.flt = new Double2DArray_2D(heightF, widthF);
			this.flt.makeDirac();

		}
	}

	public Filter_2D(ImagePlus imgF, FilterType_2D tF, int widthF, int heightF, int s) {
		this.img = IJ.createImage("Filterzzz", "32-bit", width, height, 1);
		this.img = new ImagePlus("Filter", imgF.getProcessor());
		this.width = widthF;
		this.height = heightF;
		this.type = tF;
		this.p = 2;
		this.alpha = 1;
		this.slices = s;
		if (tF == FilterType_2D.GABOR3D || s > 1 && tF == FilterType_2D.DIRAC)
			this.stackFilter = new Double3DArray_2D(s, heightF, widthF);
		else
			this.flt = new Double2DArray_2D(heightF, widthF);
		if (tF == FilterType_2D.DIRAC) {
			if (s > 1)
				this.stackFilter.makeDirac();
			else
				this.flt.makeDirac();
		}
	}

	public Filter_2D(FilterType_2D tF, int widthF, int heightF, int s) {
		this.img = IJ.createImage("Filterzzz", "32-bit", widthF, heightF, 1);
		this.width = widthF;
		this.height = heightF;
		this.setSlices(s);
		this.type = tF;
		this.p = 2;
		this.alpha = 1;
		if (tF == FilterType_2D.GABOR3D || s > 1 && tF == FilterType_2D.DIRAC) // if
																			// the
																			// filter
																			// is
																			// gabor3D
																			// OR
																			// dirac
																			// working
																			// on
																			// a
																			// stack
		{
			this.stackFilter = new Double3DArray_2D(s, heightF, widthF);
		} else {
			this.flt = new Double2DArray_2D(heightF, widthF);

		}
		if (tF == FilterType_2D.DIRAC) {
			if (s > 1) {
				this.stackFilter.makeDirac();
			} else {
				this.flt.makeDirac();
			}
		}
	}

	public Filter_2D(ImagePlus imgF, FilterType_2D typeF, int pF, double alphaF,
			int widthF, int heightF, int s) {
		this(imgF, typeF, widthF, heightF);
		this.p = pF;
		this.alpha = alphaF;
		this.setSlices(s);
		if (typeF == FilterType_2D.GABOR3D || s > 1 && typeF == FilterType_2D.DIRAC)
			this.stackFilter = new Double3DArray_2D(s, heightF, widthF);
		else
			this.flt = new Double2DArray_2D(heightF, widthF);
		if (typeF == FilterType_2D.DIRAC) {
			if (s > 1)
				this.stackFilter.makeDirac();
			else
				this.flt.makeDirac();
		}

	}

	public FilterType_2D getType() {
		return type;
	}

	public void setType(FilterType_2D type) {
		this.type = type;
	}

	public ImagePlus getImg() {
		return img;
	}

	public void setImg(ImagePlus imgF) {
		this.img = imgF;
	}

	public int getP() {
		return p;
	}

	public void setP(int p) {
		this.p = p;
	}

	public double getAlpha() {
		return alpha;
	}

	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}

	public int getWidth() {
		return width;
	}

	public void setWidth(int width) {
		this.width = width;
	}

	public int getHeight() {
		return height;
	}

	public void setHeight(int height) {
		this.height = height;
	}

	public int getSlices() {
		return slices;
	}

	public void setSlices(int slices) {
		this.slices = slices;
	}

	public Double3DArray_2D getStackFilter() {
		return stackFilter;
	}

	public void setStackFilter(Double3DArray_2D stackFilter) {
		this.stackFilter = stackFilter;
	}

	public ImageStack getStack() {
		return stack;
	}

	public void setStack(ImageStack stack) {
		this.stack = stack;
	}

	public Double2DArray_2D getFlt() {
		return flt;
	}

	public void setFlt(Double2DArray_2D flt) {
		this.flt = flt;
	}

}