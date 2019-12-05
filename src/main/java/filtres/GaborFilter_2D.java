package filtres;

import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import java.io.Serializable;
import process.Double2DArray_2D;

public class GaborFilter_2D extends Filter_2D implements Serializable {
	private static final long serialVersionUID = -7402740491860517689L;
	private double lambda;
	private double psi;
	private double gammax;
	private double gammay;
	private double theta;

	public GaborFilter_2D(ImagePlus imgF, FilterType_2D tF, int widthF, int heightF,
			int s) {
		super(imgF, tF, widthF, heightF, s);
		this.theta = 0.0D;
		this.lambda = 0.0D;
		this.psi = 0.0D;
		this.gammax = 3.0D;
		this.gammay = 1.0D;

		setAlpha(1.0D);
	}

	public GaborFilter_2D(ImagePlus imgF, FilterType_2D typeF, int pF, double alphaF,
			int widthF, int heightF, int s) {
		super(imgF, typeF, pF, alphaF, widthF, heightF, s);
		this.theta = 0.0D;
		this.lambda = 0.0D;
		this.psi = 0.0D;
		this.gammax = 3.0D;
		this.gammay = 1.0D;
		setAlpha(1.0D);
	}

	public GaborFilter_2D(ImagePlus imgF, FilterType_2D tF, int widthF, int heightF,
			double sigma, double theta, double lambda, double psi,
			double gammax, double gammay) {
		super(imgF, tF, widthF, heightF);
		this.theta = theta;
		this.lambda = lambda;
		this.psi = psi;
		this.gammax = gammax;
		this.gammay = gammay;

		setAlpha(1.0D);
	}

	public double getTheta() {
		return this.theta;
	}

	public void setTheta(double theta) {
		this.theta = theta;
	}

	public double getLambda() {
		return this.lambda;
	}

	public void setLambda(double lambda) {
		this.lambda = lambda;
	}

	public double getPsi() {
		return this.psi;
	}

	public void setPsi(double psi) {
		this.psi = psi;
	}

	public double getGamma_x() {
		return this.gammax;
	}

	public void setGamma_x(double gx) {
		this.gammax = gx;
	}

	public double getGamma_y() {
		return this.gammay;
	}

	public void setGamma_y(double gy) {
		this.gammay = gy;
	}

	public void computeNewGabor() {
		int width = getWidth();
		int height = getHeight();
		double theta = getTheta() * 3.141592653589793D / 180.0D;
		double lambda = getLambda();
		double psi = getPsi() * 3.141592653589793D / 180.0D;
		double sigmax = getGamma_x();
		double sigmay = getGamma_y();

		double alpha = Math.floor(width / 2) + 1.0D;
		double beta = Math.floor(height / 2) + 1.0D;

		double[] imgGaborD = new double[width * height];

		Double2DArray_2D temp = new Double2DArray_2D(height, width);
		ImagePlus imgGabor = IJ.createImage("Gabor Filterz", "32-bit",
				getWidth(), getHeight(), 1);
		FloatProcessor ip = (FloatProcessor) imgGabor.getProcessor();

		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {

				double x = alpha - i;
				double y = beta - j;
				double x_theta = x * Math.cos(theta) + y * Math.sin(theta);
				double y_theta = -x * Math.sin(theta) + y * Math.cos(theta);
				double pixel = Math.exp(-0.5D * x_theta * x_theta
						/ (sigmax * sigmax) - 0.5D * y_theta * y_theta
						/ (sigmay * sigmay))
						* Math.cos(x_theta * lambda / sigmax + psi);
				imgGaborD[(i + j * width)] = pixel;
			}
		}

		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				int v = i + j * width;

				ip.setf(i, j, (float) imgGaborD[v]);
				temp.setValue(imgGaborD[v], j, i, false);

				IJ.showProgress(i + j * width, width * height);
			}
			setImg(imgGabor);
			setFlt(temp);
		}
	}
}
