package process;

import java.text.DecimalFormat;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_3D;

public class Double3DArray_2D {

	private double[][][] array;
	private int rows, columns, slices;
	private DoubleFFT_3D fft;

	public Double3DArray_2D(int s, int r, int c) {
		this.slices = s;
		this.rows = r;
		this.columns = c;
		array = new double[s][r][c * 2];
	}

	public Double3DArray_2D(double[][][] tab) {
		this.slices = tab.length;
		this.rows = tab[0].length;
		this.columns = tab[0][0].length / 2;

		array = new double[slices][rows][columns * 2];

		for (int k = 0; k < slices; k++) {
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < columns * 2; j++) {
					array[k][i][j] = tab[k][i][j];
				}
			}
		}
		fft = new DoubleFFT_3D(tab.length, tab[0].length, tab[0][0].length / 2);

	}

	public double[][][] getArray() {
		return this.array;
	}

	public void setArray(double[][][] a) {
		for (int k = 0; k < slices; k++) {
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < columns * 2; j++) {
					this.array[k][i][j] = a[k][i][j];
				}
			}
		}
	}

	public void ReInitializeToZero() {
		for (int k = 0; k < slices; k++) {
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < columns * 2; j++) {
					this.array[k][i][j] = 0;
				}
			}
		}
	}

	public ImageStack getImageStack() {
		ImageStack img = IJ.createImage("Gabor Filterz", "8-bit black",
				this.getColumns(), this.getRows(), 1).getImageStack();
		ImagePlus tmp;
		img.deleteLastSlice();

		for (int k = 0; k < slices; k++) {
			tmp = IJ.createImage("Gabor Filterz " + k, "8-bit black",
					this.getColumns(), this.getRows(), 1);
			for (int i = 0; i < columns; i++) {
				for (int j = 0; j < rows; j++) {
					tmp.getProcessor().setf(i, j,
							(float) this.getValue(k, j, i, false));
				}
			}
			img.addSlice("slice #" + k, tmp.getProcessor());
		}
		return img;
	}

	public void setValue(double value, int slice, int row, int column,
			boolean complex) {
		if (slice < this.slices && row < this.rows && column < this.columns) {
			if (complex)
				this.array[slice][row][2 * column + 1] = value;
			else
				this.array[slice][row][2 * column] = value;
		} else {
			IJ.log("FATAL ERROR: indexoutofbounds SETTING value '" + value
					+ "' at" + "(" + slice + "," + row + "," + column + ")");
		}
	}

	public double getValue(int slice, int row, int column, boolean complex) {
		if (slice < this.slices && row < this.rows && column < this.columns) {
			if (complex)
				return this.array[slice][row][2 * column + 1];
			else
				return this.array[slice][row][2 * column];
		} else {
			IJ.log("FATAL ERROR: indexoutofbounds getting value at" + "("
					+ slice + "," + row + "," + column + ")");
			return 0;
		}
	}

	public int getRows() {
		return rows;
	}

	public int getColumns() {
		return columns;
	}

	public int getSlices() {
		return slices;
	}

	public void printArray() {
		if (this.rows < 10) {
			DecimalFormat f = new DecimalFormat();
			f.setMaximumFractionDigits(4);
			String str = "";
			for (int k = 0; k < array.length; k++) {
				str += "\nSlide #" + k;
				for (int i = 0; i < array[0].length; i++) {
					str += "\n-> ";
					for (int j = 0; j < array[0][0].length / 2; j++) {
						str += "(" + k + "," + i + "," + j + ")="
								+ f.format(array[k][i][2 * j]) + "_"
								+ f.format(array[k][i][2 * j + 1]) + "\t| ";
					}
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
				str += "\nSlide #" + k;
				for (int i = 0; i < array[0].length; i++) {
					str += "\n-> ";
					for (int j = 0; j < array[0][0].length / 2; j++) {
						str += "(" + k + "," + i + "," + j + ")="
								+ f.format(array[k][i][2 * j]) + "_"
								+ f.format(array[k][i][2 * j + 1]) + "\t| ";
					}
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

	// Double3DArray1=1 + Double3DArray2
	public void add1(Double3DArray_2D d) {
		double val;
		for (int i = 0; i < slices; i++) {
			for (int j = 0; j < rows; j++) {
				for (int k = 0; k < columns; k++) {
					val = 1 + d.getValue(i, j, k, false);

					this.setValue(val, i, j, k, false);
				}
			}
		}
	}

	public Double3DArray_2D multiply(Double3DArray_2D d) {
		Double3DArray_2D res = new Double3DArray_2D(this.slices, this.rows,
				this.columns);
		double[] val = new double[2];
		for (int i = 0; i < slices; i++) {
			for (int j = 0; j < rows; j++) {
				for (int k = 0; k < columns; k++) {
					val = multiplyComplex(this.getValue(i, j, k, false),
							this.getValue(i, j, k, true),
							d.getValue(i, j, k, false),
							d.getValue(i, j, k, true));

					res.setValue(val[0], i, j, k, false);
					res.setValue(val[1], i, j, k, true);
				}
			}
		}
		return res;
	}

	public void multiplyByComplexArray(Double3DArray_2D d) {

		double[] val = new double[2];
		for (int i = 0; i < slices; i++) {
			for (int j = 0; j < rows; j++) {
				for (int k = 0; k < columns; k++) {
					val = multiplyComplex(this.getValue(i, j, k, false),
							this.getValue(i, j, k, true),
							d.getValue(i, j, k, false),
							d.getValue(i, j, k, true));

					this.setValue(val[0], i, j, k, false);
					this.setValue(val[1], i, j, k, true);
				}
			}
		}
	}

	public Double3DArray_2D divide(Double3DArray_2D d) {
		Double3DArray_2D res = new Double3DArray_2D(this.slices, this.rows,
				this.columns);
		double[] val = new double[2];
		for (int i = 0; i < slices; i++) {
			for (int j = 0; j < rows; j++) {
				for (int k = 0; k < columns; k++) {
					val = dividesComplex(this.getValue(i, j, k, false),
							this.getValue(i, j, k, true),
							d.getValue(i, j, k, false),
							d.getValue(i, j, k, true));
					res.setValue(val[0], i, j, k, false);
					res.setValue(val[1], i, j, k, true);
				}
			}
		}
		return res;

	}

	public void divideByReal(double val) {
		for (int i = 0; i < slices; i++) {
			for (int j = 0; j < rows; j++) {
				for (int k = 0; k < columns; k++) {
					this.setValue(this.getValue(i, j, k, false) / val, i, j, k,
							false);
					this.setValue(this.getValue(i, j, k, true) / val, i, j, k,
							true);
				}
			}
		}
	}

	public void divideByRealArray(Double3DArray_2D d) {
		for (int i = 0; i < slices; i++) {
			for (int j = 0; j < rows; j++) {
				for (int k = 0; k < columns; k++) {
					this.setValue(
							this.getValue(i, j, k, false)
									/ d.getValue(i, j, k, false), i, j, k,
							false);
					this.setValue(
							this.getValue(i, j, k, true)
									/ d.getValue(i, j, k, false), i, j, k, true);
				}
			}
		}
	}

	public void initFFT() {
		if (this.fft == null) {
			this.fft = new DoubleFFT_3D(slices, rows, columns);
		}
	}

	public Double3DArray_2D getFFTn() {
		this.initFFT();
		// lets backup array
		double[][][] backup = new double[slices][rows][columns * 2];

		for (int k = 0; k < slices; k++) {
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < columns * 2; j++) {
					backup[k][i][j] = array[k][i][j];
				}
			}
		}

		fft.complexForward(this.array);

		Double3DArray_2D res = new Double3DArray_2D(this.array);
		for (int k = 0; k < slices; k++) {
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < columns * 2; j++) {
					array[k][i][j] = backup[k][i][j];
				}
			}
		}
		return res;
	}

	public Double3DArray_2D getIFFTn(boolean scaling) {
		this.initFFT();
		// lets backup array
		double[][][] backup = new double[slices][rows][columns * 2];

		for (int k = 0; k < slices; k++) {
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < columns * 2; j++) {
					backup[k][i][j] = array[k][i][j];
				}
			}
		}

		fft.complexInverse(this.array, scaling);

		Double3DArray_2D res = new Double3DArray_2D(this.array);
		// do not modify array
		for (int k = 0; k < slices; k++) {
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < columns * 2; j++) {
					array[k][i][j] = backup[k][i][j];
				}
			}
		}
		return res;
	}

	public Double3DArray_2D getArrayPowerOf2() {
		Double3DArray_2D res = new Double3DArray_2D(this.array);
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					res.setValue(power2(this.getValue(i, j, k, false)), i, j,
							k, false);
					res.setValue(power2(this.getValue(i, j, k, true)), i, j, k,
							true);
				}
			}
		}
		return res;
	}

	// // fPsi=fPsi+abs(fPsii).^2/(alphas(i)^2);
	public void compute_fPsi(Double3DArray_2D fPsii, double alpha) {
		double value;
		for (int i = 0; i < this.getSlices(); i++) {
			for (int j = 0; j < this.getRows(); j++) {
				for (int k = 0; k < this.getColumns(); k++) {
					value = this.getValue(i, j, k, false)
							+ power2(Math.sqrt(power2(fPsii.getValue(i, j, k,
									false))
									+ power2(fPsii.getValue(i, j, k, true))))
							/ power2(alpha);
					this.setValue(value, i, j, k, false);
				}
			}
		}
	}

	// //sumPsi=sumPsi+abs(fPsii).^2*d/(alphas(i));
	public void compute_sumpsi(Double3DArray_2D psii, Double3DArray_2D d, double alpha) {
		double value;
		for (int i = 0; i < this.getSlices(); i++) {
			for (int j = 0; j < this.getRows(); j++) {
				for (int k = 0; k < this.getColumns(); k++) {
					value = this.getValue(i, j, k, false)
							+ (power2(psii.getValue(i, j, k, false)) + power2(psii
									.getValue(i, j, k, true)))
							* d.getValue(i, j, k, false) / alpha;
					this.setValue(value, i, j, k, false);
				}
			}
		}
	}

	// fPsi=fPsi*alpha^2;
	public void compute_fPsi2(double alpha) {
		double value;
		double alphaPower2 = power2(alpha);
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					value = this.getValue(i, j, k, false) * alphaPower2;
					this.setValue(value, i, j, k, false);
					value = this.getValue(i, j, k, true) * alphaPower2;
					this.setValue(value, i, j, k, true);
				}
			}
		}
	}

	// d=abs(fd1).^2+abs(fd2).^2+abs(fd3).^2;
	public void computeD(Double3DArray_2D fd1, Double3DArray_2D fd2, Double3DArray_2D fd3) {
		double value;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					value = power2(Math.sqrt(power2(fd1
							.getValue(i, j, k, false))
							+ power2(fd1.getValue(i, j, k, true))))
							+ power2(Math.sqrt(power2(fd2.getValue(i, j, k,
									false))
									+ power2(Math.abs(fd2.getValue(i, j, k,
											true)))))
							+ power2(Math.sqrt(power2(fd3.getValue(i, j, k,
									false))
									+ power2(Math.abs(fd3.getValue(i, j, k,
											true)))));
					this.setValue(value, i, j, k, false);
				}
			}
		}
	}

	// fb=sqrt(fPsi).*d.*fu0; %Complex
	public void computeFb(Double3DArray_2D fPsi, Double3DArray_2D d, Double3DArray_2D fu0) {
		double value, sqrt;
		Double3DArray_2D multi_d_fu0 = d.multiply(fu0);
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {

					sqrt = mySqrt(fPsi.getValue(i, j, k, false));

					value = sqrt * multi_d_fu0.getValue(i, j, k, false);
					this.setValue(value, i, j, k, false);

					value = sqrt * multi_d_fu0.getValue(i, j, k, true);
					this.setValue(value, i, j, k, true);

				}
			}
		}

	}

	// flambda=fb./(d.*fPsi+alpha); %Complex
	public Double3DArray_2D computeFlambda(Double3DArray_2D fb, Double3DArray_2D d,
			Double3DArray_2D fPsi, double alpha) {
		Double3DArray_2D multi_d_fPsi = d.multiply(fPsi);
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					multi_d_fPsi.setValue(
							this.getValue(i, j, k, false) + alpha, i, j, k,
							false);
				}
			}
		}
		return fb.divide(multi_d_fPsi);
	}

	// psi=ifftn(sqrt(fPsi)); %Real
	public void computePsi(Double3DArray_2D fPsi) {
		this.initFFT();
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					double sqrt = mySqrt(fPsi.getValue(i, j, k, false));
					this.setValue(sqrt, i, j, k, false);

				}
			}
		}
		fft.complexInverse(this.array, true);
	}

	public void makeDirac() {
		this.setValue(1, 0, 0, 0, false);
	}

	// b=ifftn(flambda.*sqrt(fPsi)); %Real
	public void computeBIFFT(Double3DArray_2D fPsi, Double3DArray_2D flambda) {
		double value;
		this.initFFT();
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					double sqrt = mySqrt(fPsi.getValue(i, j, k, false));
					value = flambda.getValue(i, j, k, false) * sqrt;
					this.setValue(value, i, j, k, false);
				}
			}
		}
		fft.complexInverse(this.array, true);
	}

	// fnoise=fnoise+bruiths(i));
	public void compute_fnoise(Double3DArray_2D bruith) {
		double value;
		{
			for (int i = 0; i < this.getSlices(); i++) {
				for (int j = 0; j < this.getRows(); j++) {
					for (int k = 0; k < this.getColumns(); k++) {
						value = this.getValue(i, j, k, false)
								+ bruith.getValue(i, j, k, false);
						this.setValue(value, i, j, k, false);
						value = this.getValue(i, j, k, true)
								+ bruith.getValue(i, j, k, true);
						this.setValue(value, i, j, k, true);
					}
				}
			}
		}
	}

	public void denoiseImage(Double3DArray_2D u0, Double3DArray_2D noise) {
		this.initFFT();

		double value = -1;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					value = u0.getValue(i, j, k, false)
							- noise.getValue(i, j, k, false);
					this.setValue(value, i, j, k, false);
				}

			}
		}
	}

	public double mySqrt(double real) {
		return Math.sqrt(real);
	}

	public Double3DArray_2D getAbs() {
		Double3DArray_2D res = new Double3DArray_2D(this.slices, this.rows,
				this.columns);
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					res.setValue(
							Math.sqrt(power2(this.getValue(i, j, k, false))
									+ power2(this.getValue(i, j, k, true))), i,
							j, k, false);
				}
			}
		}
		return res;
	}

	public double getAbs(int slice, int row, int column) {
		return Math.sqrt(power2(this.getValue(slice, row, column, false))
				+ power2(this.getValue(slice, row, column, true)));
	}

	public double getNorm() {
		double sumsq = 0;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					sumsq += array[i][j][k * 2] * array[i][j][k * 2]
							+ array[i][j][k * 2 + 1] * array[i][j][k * 2 + 1];
				}
			}
		}

		return Math.sqrt(sumsq);

	}

	public Double3DArray_2D getConj() {
		Double3DArray_2D res = new Double3DArray_2D(slices, rows, columns);
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					res.setValue(this.getValue(i, j, k, false), i, j, k, false);
					res.setValue(-this.getValue(i, j, k, true), i, j, k, true);
				}
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
	public void compute_hh(Double3DArray_2D fpsi) {
		double val;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < rows; j++) {
				for (int k = 0; k < columns; k++) {
					val = fpsi.getAbs(i, j, k);
					this.setValue(val * val, i, j, k, false);
				}
			}
		}
	}

	// h=max(abs(hh(:).*fd1(:)));
	public double compute_h(Double3DArray_2D d) {
		Double3DArray_2D res = new Double3DArray_2D(slices, rows, columns);
		res = this.multiply(d);
		double max = Double.NEGATIVE_INFINITY;
		double abs;

		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					abs = res.getAbs(i, j, k);
					if (max < abs) {
						max = abs;
					}
				}
			}
		}
		return max;
	}

	// alphas=sqrt(numel(u0))*NormesInf./(norm(u0(:))*etas);
	public double[] compute_alphas(double[] etas, double[] normInf) {
		double sqrt = Math.sqrt(slices * rows * columns);
		double normu0 = this.getNorm();
		double[] alphas = new double[etas.length];
		for (int i = 0; i < etas.length; i++) {
			alphas[i] = (sqrt * normInf[i]) / (normu0 * etas[i]);
		}
		return alphas;
	}

	// fpsi=fpsi+abs(fPsis(:,:,:,i)).^2/alphas(i);
	public void compute_fPsi3(Double3DArray_2D fpsis, double alpha) {
		double val;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					val = this.getValue(i, j, k, false)
							+ power2(fpsis.getAbs(i, j, k)) / alpha;
					this.setValue(val, i, j, k, false);
				}
			}
		}
	}

	// fpsi=sqrt(alpha*fpsi);
	public void compute_fPsi4(int alpha) {
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					this.setValue(
							Math.sqrt(this.getValue(i, j, k, false) * alpha),
							i, j, k, false);

					this.setValue(
							Math.sqrt(this.getValue(i, j, k, true) * alpha), i,
							j, k, true);
				}
			}
		}
	}

	// d1u0=ifftn(fd1.*fu0);
	public Double3DArray_2D compute_dXu0(Double3DArray_2D d, Double3DArray_2D fu0) {
		Double3DArray_2D res = new Double3DArray_2D(slices, rows, columns);
		res = fu0.multiply(d);
		return res.getIFFTn(true);
	}

	// L=2/alpha*max(abs(fpsi(:)).^2.*(abs(fd1(:)).^2+abs(fd2(:)).^2+abs(fd3(:)).^2));
	public double compute_L(int alpha, Double3DArray_2D d1, Double3DArray_2D d2,
			Double3DArray_2D d3) {
		// Double3DArray res = new Double3DArray(slices, rows, columns);
		double value, max = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					value = power2(this.getAbs(i, j, k))
							* (power2(d1.getAbs(i, j, k))
									+ power2(d2.getAbs(i, j, k)) + power2(d3
										.getAbs(i, j, k)));
					if (max < value) {
						max = value;
					}
				}
			}
		}
		return 2 / (alpha) * max;
	}

	// fAy=conj(fpsi).*(fy1.*conj(fd1)+fy2.*conj(fd2)+fy3.*conj(fd3));
	public Double3DArray_2D compute_fAy(Double3DArray_2D fPsi, Double3DArray_2D fy1,
			Double3DArray_2D fy2, Double3DArray_2D fy3, Double3DArray_2D d1,
			Double3DArray_2D d2, Double3DArray_2D d3) {

		Double3DArray_2D mult1 = fy1.multiply(d1.getConj());
		Double3DArray_2D mult2 = fy2.multiply(d2.getConj());
		Double3DArray_2D mult3 = fy3.multiply(d3.getConj());
		Double3DArray_2D res = new Double3DArray_2D(slices, rows, columns);

		double val;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					val = mult1.getValue(i, j, k, false)
							+ mult2.getValue(i, j, k, false)
							+ mult3.getValue(i, j, k, false);
					res.setValue(val, i, j, k, false);

					val = mult1.getValue(i, j, k, true)
							+ mult2.getValue(i, j, k, true)
							+ mult3.getValue(i, j, k, true);
					res.setValue(val, i, j, k, true);
				}
			}
		}
		return res.multiply(fPsi.getConj());
	}

	// CF(i)=sum(d1u0(:).*y1(:)+d2u0(:).*y2(:)+d3u0(:).*y3(:))-1/(2*alpha)*norm(Ay(:))^2;
	public double compute_CF(Double3DArray_2D d1u0, Double3DArray_2D y1,
			Double3DArray_2D d2u0, Double3DArray_2D y2, Double3DArray_2D d3u0,
			Double3DArray_2D y3, int alpha) {
		double a = (double) alpha;
		Double3DArray_2D mult1 = d1u0.multiply(y1);
		Double3DArray_2D mult2 = d2u0.multiply(y2);
		Double3DArray_2D mult3 = d3u0.multiply(y3);
		double norm = this.getNorm();
		double aYnormSquared = power2(norm);
		double sum = 0;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					sum += mult1.getValue(i, j, k, false)
							+ mult2.getValue(i, j, k, false)
							+ mult3.getValue(i, j, k, false)
							+ mult1.getValue(i, j, k, true)
							+ mult2.getValue(i, j, k, true)
							+ mult3.getValue(i, j, k, true);
				}
			}
		}
		double val = sum - (1 / (2 * a)) * aYnormSquared;
		return val;
	}

	// nablaF1=d1u0-1/alpha*ifftn(fd1.*fpsi.*fAy);
	public Double3DArray_2D compute_nablaFX(int alpha, Double3DArray_2D d1,
			Double3DArray_2D fPsi, Double3DArray_2D fAy) {
		Double3DArray_2D mult1_ifftn = d1.multiply(fPsi.multiply(fAy)).getIFFTn(
				true);
		Double3DArray_2D res = new Double3DArray_2D(slices, rows, columns);

		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					res.setValue(this.getValue(i, j, k, false) - (1 / alpha)
							* (mult1_ifftn.getValue(i, j, k, false)), i, j, k,
							false);
					res.setValue(this.getValue(i, j, k, true) - (1 / alpha)
							* (mult1_ifftn.getValue(i, j, k, true)), i, j, k,
							true);

				}
			}
		}
		return res;
	}

	// q1=y1+1/L*nablaF1;
	public Double3DArray_2D compute_qX(Double3DArray_2D y1, double l,
			Double3DArray_2D nablaF) {
		Double3DArray_2D res = new Double3DArray_2D(slices, rows, columns);

		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					res.setValue((double) y1.getValue(i, j, k, false) + (1 / l)
							* nablaF.getValue(i, j, k, false), i, j, k, false);
					res.setValue((double) y1.getValue(i, j, k, true) + (1 / l)
							* nablaF.getValue(i, j, k, true), i, j, k, true);
				}
			}
		}
		return res;
	}

	// nq=sqrt(q1.^2+q2.^2+q3.^2);
	public Double3DArray_2D compute_nq(Double3DArray_2D q2, Double3DArray_2D q3) {
		Double3DArray_2D res = new Double3DArray_2D(slices, rows, columns);
		double val;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					val = Math.sqrt(power2(this.getValue(i, j, k, false))
							+ power2(q2.getValue(i, j, k, false))
							+ power2(q3.getValue(i, j, k, false)));
					res.setValue(val, i, j, k, false);
				}
			}
		}
		return res;
	}

	// q1(nq>1)=q1(nq>1)./nq(nq>1);
	public void compute_qX2(Double3DArray_2D nq) {

		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					if (nq.getValue(i, j, k, false) > 1) {
						this.setValue(
								this.getValue(i, j, k, false)
										/ nq.getValue(i, j, k, false), i, j, k,
								false);
					}
				}
			}
		}
	}

	// y1=q1+(i-1)/(i+2)*(q1-q1p);
	public Double3DArray_2D compute_yX(Double3DArray_2D q1, int inc, Double3DArray_2D q1p) {
		Double3DArray_2D res = new Double3DArray_2D(slices, rows, columns);
		double val;
		double incr = (double) inc;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					val = ((incr - 1) / (incr + 2) * (q1.getValue(i, j, k,
							false) - q1p.getValue(i, j, k, false)))
							+ q1.getValue(i, j, k, false);
					res.setValue(val, i, j, k, false);

					val = ((incr - 1) / (incr + 2) * (q1
							.getValue(i, j, k, true) - q1p.getValue(i, j, k,
							true)))
							+ q1.getValue(i, j, k, true);
					res.setValue(val, i, j, k, true);
				}
			}
		}
		return res;
	}

	// u=u0-ifftn(fpsi.*fAy)/alpha;
	public void compute_u(Double3DArray_2D u0, Double3DArray_2D fPsi,
			Double3DArray_2D fAy, int alpha, Double3DArray_2D noise) {
		Double3DArray_2D mult_ifftn = fPsi.multiply(fAy).getIFFTn(true);
		double val;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					val = u0.getValue(i, j, k, false)
							- mult_ifftn.getValue(i, j, k, false) / alpha;
					this.setValue(val, i, j, k, false);

					noise.setValue(mult_ifftn.getValue(i, j, k, false) / alpha,
							i, j, k, false);
				}
			}
		}
	}

	public boolean equals(Double3DArray_2D a) {
		if (this.getSlices() == a.getSlices() && this.getRows() == a.getRows()
				&& this.getColumns() == a.getColumns()) {
			for (int i = 0; i < array.length; i++) {
				for (int j = 0; j < array[0].length; j++) {
					for (int k = 0; k < array[0][0].length / 2; k++) {
						if ((this.getValue(i, j, k, false) != a.getValue(i, j,
								k, false))
								|| (this.getValue(i, j, k, true) != a.getValue(
										i, j, k, true)))
							return false;

					}
				}
			}
			return true;
		} else {
			return false;
		}
	}

	// return the minimum and maximum (not complex) value of the array. min is
	// stored in [0] and max in [1].
	public double[] minNmax() {
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		double val = 0;
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					val = this.getValue(i, j, k, false);
					if (val > max)
						max = val;
					if (val < min)
						min = val;
				}
			}
		}
		double[] ret = { min, max };
		return ret;
	}

	// initialise alpha for the DenoiseH1 process
	public double[] init_alphas(Double3DArray_2D fu0, Double3DArray_2D d,
			Double3DArray_2D[] psis) {
		double[] res;
		res = new double[array.length];
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[0].length; j++) {
				for (int k = 0; k < array[0][0].length / 2; k++) {
					res[i] = 0;
				}
			}
		}
		return res;
	}
}
