package process;

import java.text.DecimalFormat;
import ij.*;

public class Vsnr_denoiser_2D {

	/**
	 * 
	 * This function is denoising an image using the TV algorithm
	 * 
	 * @param u0
	 *            image to denoise
	 * @param psis
	 *            array of filters
	 * @param etas
	 *            array of noise level
	 * @param nit
	 *            number of iterations
	 * @param noise
	 *            the final noise array
	 * @param slice
	 *            number of the current slices of the stack
	 * @param nbslices
	 *            total number of slices
	 * @return the denoised image
	 */
	public static Double2DArray_2D denoiseTV_2D(Double2DArray_2D u0,
			Double2DArray_2D[] psis, double[] etas, int nit, Double2DArray_2D noise,
			int slice, int nbslices) {

		try {
			IJ.showStatus("Starting denoising ...");
			int alpha = 1;

			int nbRows = u0.getRows();
			int nbCols = u0.getColumns();

			Double2DArray_2D d1 = new Double2DArray_2D(nbRows, nbCols);
			d1.setValue(1, 0, 0, false);
			d1.setValue(-1, d1.getRows() - 1, 0, false);

			Double2DArray_2D d2 = new Double2DArray_2D(nbRows, nbCols);
			d2.setValue(1, 0, 0, false);
			d2.setValue(-1, 0, d2.getColumns() - 1, false);

			d1 = d1.getFFTn();
			d2 = d2.getFFTn();

			Double2DArray_2D fu0 = u0.getFFTn();

			Double2DArray_2D[] fPsis = new Double2DArray_2D[psis.length];
			for (int i = 0; i < psis.length; i++) {
				fPsis[i] = psis[i].getFFTn();
			}
			double[] normesInf = new double[psis.length];
			double h = 0;
			Double2DArray_2D hh = new Double2DArray_2D(nbRows, nbCols);
			for (int i = 0; i < psis.length; i++) {
				normesInf[i] = 0;
				hh.compute_hh(fPsis[i]);
				h = hh.compute_h(d1);
				if (normesInf[i] < h) {
					normesInf[i] = h;
				}

				h = hh.compute_h(d2);
				if (normesInf[i] < h) {
					normesInf[i] = h;
				}
			}
			double[] alphas = u0.compute_alphas(etas, normesInf);

			Double2DArray_2D fPsi = new Double2DArray_2D(nbRows, nbCols);

			for (int i = 0; i < psis.length; i++) {
				fPsi.compute_fPsi3(fPsis[i], alphas[i]);
			}

			fPsi.compute_fPsi4(alpha);
			double[] CF = new double[nit];
			Double2DArray_2D y1 = new Double2DArray_2D(nbRows, nbCols);
			Double2DArray_2D y2 = new Double2DArray_2D(nbRows, nbCols);

			Double2DArray_2D q1 = new Double2DArray_2D(nbRows, nbCols);
			Double2DArray_2D q2 = new Double2DArray_2D(nbRows, nbCols);

			Double2DArray_2D d1u0 = new Double2DArray_2D(nbRows, nbCols);
			d1u0 = d1u0.compute_dXu0(d1, fu0);

			Double2DArray_2D d2u0 = new Double2DArray_2D(nbRows, nbCols);
			d2u0 = d2u0.compute_dXu0(d2, fu0);
			double L = -1;

			L = fPsi.compute_L(alpha, d1, d2);
			Double2DArray_2D fy1, fy2, fAy = new Double2DArray_2D(nbRows, nbCols), aY = new Double2DArray_2D(
					nbRows, nbCols), nablaF1, nablaF2, q1p, q2p, nq = new Double2DArray_2D(
					nbRows, nbCols);

			DecimalFormat f = new DecimalFormat();
			f.setMaximumFractionDigits(2);
			for (int i = 1; i <= nit; i++) {

				try {
					fy1 = y1.getFFTn();
					fy2 = y2.getFFTn();
					fAy = fAy.compute_fAy(fPsi, fy1, fy2, d1, d2);
					aY = fAy.getIFFTn(true);

					CF[i - 1] = aY.compute_CF(d1u0, y1, d2u0, y2, alpha);

					nablaF1 = d1u0.compute_nablaFX(alpha, d1, fPsi, fAy);
					nablaF2 = d2u0.compute_nablaFX(alpha, d2, fPsi, fAy);
					q1p = new Double2DArray_2D(q1.getArray());
					q2p = new Double2DArray_2D(q2.getArray());

					q1 = q1.compute_qX(y1, L, nablaF1);
					q2 = q2.compute_qX(y2, L, nablaF2);

					nq = q1.compute_nq(q2);

					q1.compute_qX2(nq);
					q2.compute_qX2(nq);

					y1 = y1.compute_yX(q1, i, q1p);
					y2 = y2.compute_yX(q2, i, q2p);

					System.gc();
				} catch (NullPointerException npe) {
					y1 = null;
					y2 = null;
					d1 = null;
					d2 = null;
					nq = null;
					q1 = null;
					q2 = null;
					q1p = null;
					q2p = null;
					nablaF1 = null;
					nablaF2 = null;
					d1u0 = null;
					d2u0 = null;
					fy1 = null;
					fy2 = null;
					fAy = null;
					aY = null;
					hh = null;
					CF = null;
					fu0 = null;
					fPsi = null;
					fPsis = null;
					System.gc();
				}
				// we warn the user that we are denoising
				IJ.showStatus("Denoising... Iteration n" + i);
				IJ.showProgress(i, (nit + 1));
			}

			Double2DArray_2D u = new Double2DArray_2D(nbRows, nbCols);

			u.compute_u(u0, fPsi, fAy, alpha, noise);

			// We should clear the memory
			y1 = null;
			y2 = null;
			d1 = null;
			d2 = null;
			nq = null;
			q1 = null;
			q2 = null;
			q1p = null;
			q2p = null;
			nablaF1 = null;
			nablaF2 = null;
			d1u0 = null;
			d2u0 = null;
			fy1 = null;
			fy2 = null;
			fAy = null;
			aY = null;
			hh = null;
			CF = null;
			fu0 = null;
			fPsi = null;
			fPsis = null;
			System.gc();

			// everything is fine, so we return the denoised image !
			IJ.showStatus("Denoising ended !");
			return u;
		} catch (Exception e) {
			System.gc();
			IJ.handleException(e);
			return null;
		}
	}

	/**
	 * 
	 * This function is denoising an image using the TV algorithm
	 * 
	 * @param u0
	 *            image to denoise
	 * @param psis
	 *            array of filters
	 * @param etas
	 *            array of noise level
	 * @param noise
	 *            the final noise array
	 * @return
	 */
	public static Double2DArray_2D denoiseH1_2D(Double2DArray_2D u0,
			Double2DArray_2D[] psis, double[] etas, Double2DArray_2D noise) {
		try {
			// we warn the user that we are denoising
			IJ.showStatus("Starting denoising ...");
			int nbRows = u0.getRows();
			int nbCols = u0.getColumns();

			Double2DArray_2D d1 = new Double2DArray_2D(nbRows, nbCols);
			d1.setValue(1, 0, 0, false);
			d1.setValue(-1, d1.getRows() - 1, 0, false);

			Double2DArray_2D d2 = new Double2DArray_2D(nbRows, nbCols);
			d2.setValue(1, 0, 0, false);
			d2.setValue(-1, 0, d2.getColumns() - 1, false);

			d1 = d1.getFFTn();
			d2 = d2.getFFTn();

			Double2DArray_2D fu0 = new Double2DArray_2D(nbRows, nbCols);
			fu0 = u0.getFFTn(); // fu0=fftn(u0);
			double norm_fu0 = fu0.getNorm();

			Double2DArray_2D[] fPsis = new Double2DArray_2D[psis.length];
			for (int i = 0; i < psis.length; i++) {
				fPsis[i] = psis[i].getFFTn();
			}

			int m = psis.length;
			double[] nbruits;
			nbruits = new double[psis.length]; // morozov's estimated noises

			for (int i = 0; i < psis.length; i++) {

				nbruits[i] = etas[i] * norm_fu0;// each component of the Morozov
												// estimated noises

			}

			Double2DArray_2D d = new Double2DArray_2D(nbRows, nbCols);
			d.computeD(d1, d2);
			d1 = null;
			d2 = null;

			Double2DArray_2D temp = new Double2DArray_2D(nbRows, nbCols);
			double[] alphas = new double[m];

			for (int i = 0; i < m; i++) {
				temp = (fPsis[i].getAbs()).getArrayPowerOf2();
				temp.multiplyByComplexArray(d);
				temp.multiplyByComplexArray(fu0);
				alphas[i] = temp.getNorm();
				alphas[i] = alphas[i] / nbruits[i];
			}

			Double2DArray_2D sumpsi = new Double2DArray_2D(nbRows, nbCols);
			Double2DArray_2D[] bruiths = new Double2DArray_2D[m];
			double[] nbruiths = new double[m];
			for (int i = 0; i < m; i++) {
				sumpsi.compute_sumpsi(fPsis[i], d, alphas[i]);
			}

			Double2DArray_2D d_multiply_fu0 = d.multiply(fu0);
			temp.add1(sumpsi);
			for (int i = 0; i < m; i++) {
				bruiths[i] = new Double2DArray_2D(nbRows, nbCols);
				bruiths[i] = (fPsis[i].getAbs()).getArrayPowerOf2();
				bruiths[i].multiplyByComplexArray(d_multiply_fu0);
				bruiths[i].divideByRealArray(temp);
				nbruiths[i] = (bruiths[i].getNorm()) / alphas[i];
			}

			boolean b = true;
			int nit = 0;
			int maxiter = 5;
			while (b && nit <= maxiter) {
				alphas[0] = alphas[0] * nbruiths[0] / nbruits[0];
				for (int i = 0; i < m - 1; i++) {
					sumpsi.ReInitializeToZero();
					for (int j = 0; j < m; j++) {
						sumpsi.compute_sumpsi(fPsis[j], d, alphas[j]);
					}
					bruiths[i + 1] = (fPsis[i + 1].getAbs()).getArrayPowerOf2();
					bruiths[i + 1].multiplyByComplexArray(d_multiply_fu0);
					temp.add1(sumpsi);
					bruiths[i + 1].divideByRealArray(temp);
					nbruiths[i + 1] = (bruiths[i + 1].getNorm())
							/ alphas[i + 1];
					alphas[i + 1] = (nbruiths[i + 1] / nbruits[i + 1])
							* alphas[i + 1];
				}
				sumpsi.ReInitializeToZero();
				for (int j = 0; j < m; j++) {
					sumpsi.compute_sumpsi(fPsis[j], d, alphas[j]);
				}

				temp.add1(sumpsi);
				for (int i = 0; i < m; i++) {
					bruiths[i] = (fPsis[i].getAbs()).getArrayPowerOf2();
					bruiths[i].multiplyByComplexArray(d_multiply_fu0);
					bruiths[i].divideByReal(alphas[i]);
					bruiths[i].divideByRealArray(temp);
					nbruiths[i] = bruiths[i].getNorm();
				}
				int i = 0;
				int j = 0;
				while (i < m && b) {
					if ((nbruits[i] - nbruiths[i]) / nbruits[i] < (1 / 10)) {
						j++;
					}
					i++;
				}
				if (j == m)
					b = false;

				// we warn the user that we are denoising
				IJ.showStatus("Denoising ...");
				IJ.showProgress(nit, (maxiter + 1));
				nit++;

			}
			Double2DArray_2D fnoise = new Double2DArray_2D(nbRows, nbCols);
			for (int i = 0; i < m; i++) {
				fnoise.compute_fnoise(bruiths[i]);
			}
			fnoise = fnoise.getIFFTn(true);
			noise.setArray(fnoise.getArray());

			Double2DArray_2D u = new Double2DArray_2D(nbRows, nbCols);
			u.denoiseImage(u0, noise);

			IJ.showStatus("Denoising done !");
			return u;
		} catch (Exception e) {
			IJ.handleException(e);
		}
		return null;
	}
}
