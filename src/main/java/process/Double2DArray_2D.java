package process;

import java.text.DecimalFormat;
import ij.IJ;
import ij.ImagePlus;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_2D;


public class Double2DArray_2D {

	private double[][] array;
	private int rows, columns;
	private DoubleFFT_2D fft;

	public Double2DArray_2D(int r, int c) {
		this.rows = r;
		this.columns = c;
		try {
			array = new double[r][c * 2];
			for (int i = 0; i < r; i++) {
				for (int j = 0; j < 2 * c; j++) {
					this.array[i][j] = 0;
				}

			}

		} catch (Throwable e) {
			IJ.error("Fiji's running OUT OF MEMORY !\n You'd better restart it NOW... !");
		}
	}

	public Double2DArray_2D(double[][] tab) {
		this.rows = tab.length;
		this.columns = tab[0].length / 2;
		boolean stop = false;
		try {
			array = new double[rows][columns * 2];

		} catch (Throwable e) {
			IJ.error("Fiji's running OUT OF MEMORY !\n You'd better restart it NOW... !");

		}
		try {
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < columns * 2; j++) {
					array[i][j] = tab[i][j];
				}
			}
			fft = new DoubleFFT_2D(tab.length, tab[0].length / 2);
		} catch (NullPointerException npe) {
			if (stop)
				return;
			else
				IJ.handleException(npe);
		}
	}

	public void ReInitializeToZero() {

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				this.array[i][2 * j] = 0;
				this.array[i][2 * j + 1] = 0;
			}
		}

	}

	public double[][] getArray() {
		return this.array;
	}

	public void setArray(double[][] a) {
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns * 2; j++) {
				this.array[i][j] = a[i][j];
			}
		}
	}

	public ImagePlus getImagePlus() {
		ImagePlus tmp = IJ.createImage("Image 2D", "8-bit black",
				this.getColumns(), this.getRows(), 1);
		for (int i = 0; i < columns; i++) {
			for (int j = 0; j < rows; j++) {
				tmp.getProcessor().setf(i, j,
						(float) this.getValue(j, i, false));
			}
		}
		return tmp;
	}

	public void setValue(double value, int row, int column, boolean complex) {
		if (row < this.rows && column < this.columns) {
			if (complex)
				this.array[row][2 * column + 1] = value;
			else
				this.array[row][2 * column] = value;
		} else {
			IJ.log("FATAL ERROR: indexoutofbounds SETTING value '" + value
					+ "' at" + "(" + row + "," + column + ")");
			try {
				throw new Exception("merdouille");
			} catch (Exception e) {
				IJ.handleException(e);
			}
		}
	}

	public double getValue(int row, int column, boolean complex) {
		if (row < this.rows && column < this.columns) {
			if (complex)
				return this.array[row][2 * column + 1];
			else
				return this.array[row][2 * column];
		} else {
			IJ.log("FATAL ERROR: indexoutofbounds getting value at" + "(" + row
					+ "," + column + ")");
			return 0;
		}
	}

	public int getRows() {
		return rows;
	}

	public int getColumns() {
		return columns;
	}

	public void printArray() {
		if (this.rows < 10) {
			DecimalFormat f = new DecimalFormat();
			f.setMaximumFractionDigits(4);
			String str = "";
			for (int k = 0; k < array.length; k++) {
				for (int i = 0; i < array[0].length / 2; i++) {
					str += "(" + k + "," + i + ")=" + f.format(array[k][i * 2])
							+ "_" + f.format(array[k][i * 2]) + "\t| ";
				}
			}
			IJ.log(str);
		}
	}

	public void printArray(String str1) {
		if (this.rows < 10) {
			DecimalFormat f = new DecimalFormat();
			f.setMaximumFractionDigits(4);
			String str = "\nPrinting _" + str1 + " :";
			for (int k = 0; k < array.length; k++) {
				for (int i = 0; i < array[0].length / 2; i++) {
					str += "(" + k + "," + i + ")=" + f.format(array[k][i * 2])
							+ "_" + f.format(array[k][i * 2 + 1]) + "\t| ";
				}
			}
			IJ.log(str + "\n");
		}
	}

	public void printArray(double[][][] arr) {
		if (arr[0][0].length < 10) {
			DecimalFormat f = new DecimalFormat();
			f.setMaximumFractionDigits(2);
			String str = "";
			for (int k = 0; k < arr.length; k++) {
				str += "\nSlide #" + k;
				for (int i = 0; i < arr[0].length; i++) {
					str += "\n-> ";
					for (int j = 0; j < arr[0][0].length / 2; j++) {
						str += "(" + k + "," + i + "," + j + ")="
								+ f.format(arr[k][i][2 * j]) + "_"
								+ f.format(arr[k][i][2 * j + 1]) + "\t| ";
					}
				}
			}
			IJ.log(str);
		}
	}

	public static double[] multiplyComplex(double real1, double im1,
			double real2, double im2) {
		double r = real1 * real2 - im1 * im2;
		double imag = real1 * im2 + im1 * real2;
		double[] res = { r, imag };
		return res;
	}

	public static double[] reciprocalComplex(double re, double im) {
		double scale = re * re + im * im;
		double[] res = { re / scale, -im / scale };
		return res;
	}

	public static double[] dividesComplex(double re1, double im1, double re2,
			double im2) {
		double[] res = reciprocalComplex(re2, im2);
		return multiplyComplex(re1, im1, res[0], res[1]);
	}

	public Double2DArray_2D multiply(Double2DArray_2D d) {
		Double2DArray_2D res = new Double2DArray_2D(this.rows, this.columns);
		double[] val = new double[2];
		for (int j = 0; j < rows; j++) {
			for (int k = 0; k < columns; k++) {
				val = multiplyComplex(this.getValue(j, k, false),
						this.getValue(j, k, true), d.getValue(j, k, false),
						d.getValue(j, k, true));

				res.setValue(val[0], j, k, false);
				res.setValue(val[1], j, k, true);
			}
		}
		return res;
	}

	public void multiplyByComplexArray(Double2DArray_2D d) {
		double[] val = new double[2];
		for (int j = 0; j < rows; j++) {
			for (int k = 0; k < columns; k++) {
				val = multiplyComplex(this.getValue(j, k, false),
						this.getValue(j, k, true), d.getValue(j, k, false),
						d.getValue(j, k, true));

				this.setValue(val[0], j, k, false);
				this.setValue(val[1], j, k, true);
			}
		}
	}

	public Double2DArray_2D divide(Double2DArray_2D d) {
		Double2DArray_2D res = new Double2DArray_2D(this.rows, this.columns);
		double[] val = new double[2];
		for (int j = 0; j < rows; j++) {
			for (int k = 0; k < columns; k++) {
				val = dividesComplex(this.getValue(j, k, false),
						this.getValue(j, k, true), d.getValue(j, k, false),
						d.getValue(j, k, true));
				res.setValue(val[0], j, k, false);
				res.setValue(val[1], j, k, true);
			}
		}
		return res;
	}

	public void divideByReal(double val) {
		for (int j = 0; j < rows; j++) {
			for (int k = 0; k < columns; k++) {
				this.setValue(this.getValue(j, k, true) / val, j, k, true);
				this.setValue(this.getValue(j, k, false) / val, j, k, false);
			}
		}
	}

	public void divideByRealArray(Double2DArray_2D d) {
		for (int j = 0; j < rows; j++) {
			for (int k = 0; k < columns; k++) {
				this.setValue(
						this.getValue(j, k, false) / d.getValue(j, k, false),
						j, k, false);
				this.setValue(
						this.getValue(j, k, true) / d.getValue(j, k, false), j,
						k, true);
			}
		}
	}

	public void initFFT() {
		if (this.fft == null) {
			this.fft = new DoubleFFT_2D(rows, columns);
		}
	}

	public Double2DArray_2D getFFTn() {
		try {
			this.initFFT();
			// lets backup array
			double[][] backup = new double[rows][columns * 2];

			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < columns * 2; j++) {
					backup[i][j] = array[i][j];
				}
			}

			fft.complexForward(this.array);

			Double2DArray_2D res = new Double2DArray_2D(this.array);
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < columns * 2; j++) {
					array[i][j] = backup[i][j];
				}
			}
			return res;
		} catch (Exception e) {
			IJ.handleException(e);
		}
		return null;
	}

	public Double2DArray_2D getIFFTn(boolean scaling) {
		this.initFFT();
		// lets backup array
		double[][] backup = new double[rows][columns * 2];

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns * 2; j++) {
				backup[i][j] = this.array[i][j];
			}
		}

		fft.complexInverse(this.array, scaling);

		Double2DArray_2D res = new Double2DArray_2D(this.array);
		// do not modify array
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns * 2; j++) {
				array[i][j] = backup[i][j];
			}
		}
		return res;
	}

	public Double2DArray_2D getArrayPowerOf2() {
		Double2DArray_2D res = new Double2DArray_2D(this.array);
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length / 2; j++) {
				res.setValue(power2(this.getValue(i, j, false)), i, j, false);
				res.setValue(power2(this.getValue(i, j, true)), i, j, true);
			}
		}
		return res;
	}

	public void makeDirac() {
		this.setValue(1, 0, 0, false);
	}

	public double mySqrt(double real) {
		return Math.sqrt(real);
	}

	public Double2DArray_2D getAbs() {
		Double2DArray_2D res = new Double2DArray_2D(this.rows, this.columns);
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length / 2; j++) {
				res.setValue(
						Math.sqrt(power2(this.getValue(i, j, false))
								+ power2(this.getValue(i, j, true))), i, j,
						false);
			}
		}
		return res;
	}

	public double getAbs(int row, int column) {
		return Math.sqrt(power2(this.getValue(row, column, false))
				+ power2(this.getValue(row, column, true)));
	}

	public double getNorm() {
		double sumsq = 0;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length / 2; j++) {
				sumsq += array[i][j * 2] * array[i][j * 2]
						+ array[i][j * 2 + 1] * array[i][j * 2 + 1];

			}
		}

		return Math.sqrt(sumsq);

	}

	public Double2DArray_2D getConj() {
		Double2DArray_2D res = new Double2DArray_2D(rows, columns);
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length / 2; j++) {
				res.setValue(this.getValue(i, j, false), i, j, false);
				res.setValue(-this.getValue(i, j, true), i, j, true);
			}
		}
		return res;
	}

	public double power2(double val) {
		return val * val;
	}

	// -------------------------------------- Pour le denoiseTV
	// -------------------------------------------

	// hh=abs(fPsis(:,:,:,i)).^2;
	public void compute_hh(Double2DArray_2D fpsi) {
		double val;
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				val = fpsi.getAbs(i, j);
				this.setValue(val * val, i, j, false);
			}
		}
	}

	// h=max(abs(hh(:).*fd1(:)));
	public double compute_h(Double2DArray_2D d) {
		Double2DArray_2D res = new Double2DArray_2D(rows, columns);
		res = this.multiply(d);
		double max = Double.NEGATIVE_INFINITY;
		double abs;

		for (int j = 0; j < array.length; j++) {
			for (int k = 0; k < array[0].length / 2; k++) {
				abs = res.getAbs(j, k);
				if (max < abs) {
					max = abs;
				}
			}
		}
		return max;
	}

	// alphas=sqrt(numel(u0))*NormesInf./(norm(u0(:))*etas);
	public double[] compute_alphas(double[] etas, double[] normInf) {
		double sqrt = Math.sqrt(rows * columns);
		double normu0 = this.getNorm();
		double[] alphas = new double[etas.length];
		for (int i = 0; i < etas.length; i++) {
			alphas[i] = (sqrt * normInf[i]) / (normu0 * etas[i]);
		}
		return alphas;
	}

	// fpsi=fpsi+abs(fPsis(:,:,:,i)).^2/alphas(i);
	public void compute_fPsi3(Double2DArray_2D fpsis, double alpha) {
		double val;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length / 2; j++) {
				val = this.getValue(i, j, false) + power2(fpsis.getAbs(i, j))
						/ alpha;
				this.setValue(val, i, j, false);
			}
		}
	}

	// fpsi=sqrt(alpha*fpsi);
	public void compute_fPsi4(int alpha) {
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length / 2; j++) {
				this.setValue(Math.sqrt(this.getValue(i, j, false) * alpha), i,
						j, false);
				this.setValue(Math.sqrt(this.getValue(i, j, true) * alpha), i,
						j, true);
			}
		}
	}

	// d1u0=ifftn(fd1.*fu0);
	public Double2DArray_2D compute_dXu0(Double2DArray_2D d, Double2DArray_2D fu0) {
		Double2DArray_2D res = new Double2DArray_2D(rows, columns);
		res = fu0.multiply(d);
		return res.getIFFTn(true);
	}

	// d=abs(fd1).^2+abs(fd2).^2;
	public void computeD(Double2DArray_2D fd1, Double2DArray_2D fd2) {
		double value;

		for (int j = 0; j < array.length; j++) {
			for (int k = 0; k < array[0].length / 2; k++) {
				value = power2(fd1.getValue(j, k, false))
						+ power2(fd1.getValue(j, k, true))
						+ power2(fd2.getValue(j, k, false))
						+ power2(Math.abs(fd2.getValue(j, k, true)));
				this.setValue(value, j, k, false);
			}
		}

	}

	// //sumPsi=sumPsi+abs(fPsii).^2*d/(alphas(i));
	public void compute_sumpsi(Double2DArray_2D psii, Double2DArray_2D d, double alpha) {
		double value;

		for (int j = 0; j < this.getRows(); j++) {
			for (int k = 0; k < this.getColumns(); k++) {
				value = this.getValue(j, k, false)
						+ (power2(psii.getValue(j, k, false)) + power2(psii
								.getValue(j, k, true)))
						* d.getValue(j, k, false) / alpha;
				this.setValue(value, j, k, false);
			}
		}

	}

	// //fnoise=fnoise+bruiths(i));
	public void compute_fnoise(Double2DArray_2D bruith) {
		double value;
		{

			for (int j = 0; j < this.getRows(); j++) {
				for (int k = 0; k < this.getColumns(); k++) {
					value = this.getValue(j, k, false)
							+ bruith.getValue(j, k, false);
					this.setValue(value, j, k, false);
					value = this.getValue(j, k, true)
							+ bruith.getValue(j, k, true);
					this.setValue(value, j, k, true);
				}
			}

		}
	}

	public void denoiseImage(Double2DArray_2D u0, Double2DArray_2D noise) {
		this.initFFT();

		double value = -1;

		for (int j = 0; j < array.length; j++) {
			for (int k = 0; k < array[0].length / 2; k++) {
				value = u0.getValue(j, k, false) - noise.getValue(j, k, false);
				this.setValue(value, j, k, false);
			}

		}
	}

	// Double2DArray1=1 + Double2DArray2
	public void add1(Double2DArray_2D d) {
		double val;

		for (int j = 0; j < rows; j++) {
			for (int k = 0; k < columns; k++) {
				val = 1 + d.getValue(j, k, false);

				this.setValue(val, j, k, false);
			}
		}

	}

	// L=2/alpha*max(abs(fpsi(:)).^2.*(abs(fd1(:)).^2+abs(fd2(:)).^2+abs(fd3(:)).^2));
	public double compute_L(int alpha, Double2DArray_2D d1, Double2DArray_2D d2) {
		double value, max = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length / 2; j++) {
				value = power2(this.getAbs(i, j))
						* (power2(d1.getAbs(i, j)) + power2(d2.getAbs(i, j)));
				if (max < value) {
					max = value;
				}
			}
		}
		return 2 / (alpha) * max;
	}

	// fAy=conj(fpsi).*(fy1.*conj(fd1)+fy2.*conj(fd2)+fy3.*conj(fd3));
	public Double2DArray_2D compute_fAy(Double2DArray_2D fPsi, Double2DArray_2D fy1,
			Double2DArray_2D fy2, Double2DArray_2D d1, Double2DArray_2D d2) {

		Double2DArray_2D mult1 = fy1.multiply(d1.getConj());
		Double2DArray_2D mult2 = fy2.multiply(d2.getConj());
		Double2DArray_2D res = new Double2DArray_2D(rows, columns);

		double val;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length / 2; j++) {
				val = mult1.getValue(i, j, false) + mult2.getValue(i, j, false);
				res.setValue(val, i, j, false);

				val = mult1.getValue(i, j, true) + mult2.getValue(i, j, true);
				res.setValue(val, i, j, true);
			}
		}
		return res.multiply(fPsi.getConj());
	}

	// CF(i)=sum(d1u0(:).*y1(:)+d2u0(:).*y2(:)+d3u0(:).*y3(:))-1/(2*alpha)*norm(Ay(:))^2;
	public double compute_CF(Double2DArray_2D d1u0, Double2DArray_2D y1,
			Double2DArray_2D d2u0, Double2DArray_2D y2, int alpha) {
		double a = (double) alpha;
		Double2DArray_2D mult1 = d1u0.multiply(y1);
		Double2DArray_2D mult2 = d2u0.multiply(y2);
		double norm = this.getNorm();
		double aYnormSquared = power2(norm);
		double sum = 0;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length / 2; j++) {
				sum += mult1.getValue(i, j, false)
						+ mult2.getValue(i, j, false)
						+ mult1.getValue(i, j, true)
						+ mult2.getValue(i, j, true);
			}
		}
		double val = sum - (1 / (2 * a)) * aYnormSquared;
		return val;
	}

	// nablaF1=d1u0-1/alpha*ifftn(fd1.*fpsi.*fAy);
	public Double2DArray_2D compute_nablaFX(int alpha, Double2DArray_2D d1,
			Double2DArray_2D fPsi, Double2DArray_2D fAy) {
		Double2DArray_2D mult1_ifftn = d1.multiply(fPsi.multiply(fAy)).getIFFTn(
				true);
		Double2DArray_2D res = new Double2DArray_2D(rows, columns);

		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length / 2; j++) {
				res.setValue(this.getValue(i, j, false) - (1 / alpha)
						* (mult1_ifftn.getValue(i, j, false)), i, j, false);
				res.setValue(this.getValue(i, j, true) - (1 / alpha)
						* (mult1_ifftn.getValue(i, j, true)), i, j, true);

			}
		}
		return res;
	}

	// q1=y1+1/L*nablaF1;
	public Double2DArray_2D compute_qX(Double2DArray_2D y1, double l,
			Double2DArray_2D nablaF) {
		Double2DArray_2D res = new Double2DArray_2D(rows, columns);

		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length / 2; j++) {
				res.setValue((double) y1.getValue(i, j, false) + (1 / l)
						* nablaF.getValue(i, j, false), i, j, false);
				res.setValue((double) y1.getValue(i, j, true) + (1 / l)
						* nablaF.getValue(i, j, true), i, j, true);
			}
		}
		return res;
	}

	// nq=sqrt(q1.^2+q2.^2);
	public Double2DArray_2D compute_nq(Double2DArray_2D q2) {
		Double2DArray_2D res = new Double2DArray_2D(rows, columns);
		double val;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length / 2; j++) {
				val = Math.sqrt(power2(this.getValue(i, j, false))
						+ power2(q2.getValue(i, j, false)));
				res.setValue(val, i, j, false);
			}
		}
		return res;
	}

	// q1(nq>1)=q1(nq>1)./nq(nq>1);
	public void compute_qX2(Double2DArray_2D nq) {

		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length / 2; j++) {
				if (nq.getValue(i, j, false) > 1) {
					this.setValue(
							this.getValue(i, j, false)
									/ nq.getValue(i, j, false), i, j, false);
				}
			}
		}
	}

	// y1=q1+(i-1)/(i+2)*(q1-q1p);
	public Double2DArray_2D compute_yX(Double2DArray_2D q1, int inc, Double2DArray_2D q1p) {
		Double2DArray_2D res = new Double2DArray_2D(rows, columns);
		double val;
		double incr = (double) inc;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length / 2; j++) {
				val = ((incr - 1) / (incr + 2) * (q1.getValue(i, j, false) - q1p
						.getValue(i, j, false))) + q1.getValue(i, j, false);
				res.setValue(val, i, j, false);

				val = ((incr - 1) / (incr + 2) * (q1.getValue(i, j, true) - q1p
						.getValue(i, j, true))) + q1.getValue(i, j, true);
				res.setValue(val, i, j, true);
			}
		}
		return res;
	}

	// u=u0-ifftn(fpsi.*fAy)/alpha;
	public void compute_u(Double2DArray_2D u0, Double2DArray_2D fPsi,
			Double2DArray_2D fAy, int alpha, Double2DArray_2D noise) {
		Double2DArray_2D mult_ifftn = fPsi.multiply(fAy).getIFFTn(true);
		double val;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length / 2; j++) {
				val = u0.getValue(i, j, false)
						- mult_ifftn.getValue(i, j, false) / alpha;
				this.setValue(val, i, j, false);

				noise.setValue(mult_ifftn.getValue(i, j, false) / alpha, i, j,
						false);
			}
		}
	}

	public boolean equals(Double2DArray_2D a) {
		if (this.getRows() == a.getRows()
				&& this.getColumns() == a.getColumns()) {
			for (int i = 0; i < array.length; i++) {
				for (int j = 0; j < array[0].length / 2; j++) {
					if ((this.getValue(i, j, false) != a.getValue(i, j, false))
							|| (this.getValue(i, j, true) != a.getValue(i, j,
									true)))
						return false;

				}
			}
			return true;
		} else {
			return false;
		}
	}

}
