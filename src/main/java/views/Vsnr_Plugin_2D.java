package views;

import filtres.FilterType_2D;
import filtres.Filter_2D;
import filtres.GaborFilter_2D;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import ij.io.OpenDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.awt.Font;
import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;

import process.Double2DArray_2D;
import process.Double3DArray_2D;
import process.Vsnr_denoiser_2D;

/**
 * 
 * @author benjamin.font
 *	
 * Plugin developer : Benjamin Font, Leo Mouly, Pierre Weiss
 * Algorithm : Jerome Fehrenbach, Pierre Weiss, Corinne Lorenzo 
 *
 *
 *         This plugin is used for denoise 2D images. To denoise stacks, uses the
 *         VSNR_3D plugin.
 */
public class Vsnr_Plugin_2D implements PlugInFilter {

	// Used for the GUI and the filters creation
	private String inputMethod;
	private String filterType;
	private ImagePlus image;
	private ImagePlus noiseImage;
	private int width;
	private int height;
	private String type;
	private ArrayList<Filter_2D> filterList = new ArrayList<Filter_2D>();

	// Filters parameters
	private double noiseLevel = 1;
	private double gammaX = 3;
	private double gammaY = 1;
	private double angle = 0;
	private double phase_psi = 0;
	private double frequency = 0;

	// Denoising parameters
	private String algorithm;
	private double nit = 50;

	// error for the text file
	private boolean errorInTxtFile = false;

	/**
	 * Run function
	 */
	@Override
	public void run(ImageProcessor ip) {
		this.image.getType();
		if (ip.equals(null)) {
			GenericDialog g = new GenericDialog("Open an image please !");
			g.setLocationRelativeTo(null);
			g.pack();
			g.setVisible(true);
		}
		if (this.image.getImageStack().getSize() > 1) {
			GenericDialog g = new GenericDialog("Error");
			g.addMessage("Use VSNR_Macro_3D for stacks !");
			g.setLocationRelativeTo(null);
			g.pack();
			g.setVisible(true);
		} else {
			if (this.askMethod()) {
				if (this.inputMethod.equals("Graphic User Interface (GUI)")) {
					if (this.askFilter()) {
						if (this.askParameters()) {
							// the user want to add more than 1 filter
							while (this.askOtherFilter()) {
								if (this.askFilter()) {
									this.askParameters();
								}
							}
							if (this.askIterationsAndAlgorithm()) {
								printParameters();
								// the algorithm is TV
								if (this.algorithm.equals("TV")) {
									// if the image is a simple 2D image
									ImageWindow wimg = new ImageWindow(
											this.image);
									wimg.setTitle("Image");
									wimg.setVisible(true);

									ImageWindow wtv2d = new ImageWindow(
											denoiseTV_2D());
									wtv2d.setTitle("Denoised Image - TV");
									wtv2d.setVisible(true);

									ImageWindow wflt = new ImageWindow(
											this.noiseImage);
									wflt.setTitle("Filtre");
									wflt.setVisible(true);
								}
								// the algorithm is H1
								else {
									// if the image is a simple 2D image
									ImageWindow wimg = new ImageWindow(
											this.image);
									wimg.setTitle("Image");
									wimg.setVisible(true);

									ImageWindow wtv2d = new ImageWindow(
											denoiseH1_2D());
									wtv2d.setTitle("Denoised Image - H1");
									wtv2d.setVisible(true);

									ImageWindow wflt = new ImageWindow(
											this.noiseImage);
									wflt.setTitle("Filtre");
									wflt.setVisible(true);
								}
							}
						}
					}
				}
				// the user want to use texts files (.txt)
				else {
					// we read the file
					if (this.readFile()) {

						// the algorithm is TV
						if (this.algorithm.equals("TV")) {
							// if the image is a simple 2D image
							ImageWindow wimg = new ImageWindow(this.image);
							wimg.setTitle("Image");
							wimg.setVisible(true);

							ImageWindow wtv2d = new ImageWindow(denoiseTV_2D());
							wtv2d.setTitle("Denoised Image - TV");
							wtv2d.setVisible(true);

							ImageWindow wflt = new ImageWindow(this.noiseImage);
							wflt.setTitle("Filtre");
							wflt.setVisible(true);
						}
						// the algorithm is H1
						else {
							// if the image is a simple 2D image
							ImageWindow wimg = new ImageWindow(this.image);
							wimg.setTitle("Image");
							wimg.setVisible(true);

							ImageWindow wtv2d = new ImageWindow(denoiseH1_2D());
							wtv2d.setTitle("Denoised Image - H1");
							wtv2d.setVisible(true);

							ImageWindow wflt = new ImageWindow(this.noiseImage);
							wflt.setTitle("Filtre");
							wflt.setVisible(true);
						}
					}
				}
			}
		}
	}

	/**
	 * setup function
	 */
	public int setup(String arg0, ImagePlus img) {
		if (img == null) {
			IJ.error("Open an image please !");
		}
		this.image = img;
		this.width = img.getWidth();
		this.height = img.getHeight();
		this.type = this.image.getBitDepth()+"-bit";
		return DOES_ALL;
	}

	/**
	 * this method print the text file parameters in a Fiji log
	 */
	public void printParameters() {
		IJ.log("Denoising_Algorithm: " + this.algorithm);
		IJ.log("Iteration_Number: " + this.nit);
		IJ.log("***");
		for (Filter_2D f : this.filterList) {
			if (f instanceof GaborFilter_2D) {
				GaborFilter_2D f2 = (GaborFilter_2D) f;
				IJ.log("Filter_Type: Gabor");
				IJ.log("Noise_Level: " + f2.getAlpha());
				IJ.log("Sigma_X: " + f2.getGamma_x());
				IJ.log("Sigma_Y: " + f2.getGamma_y());
				IJ.log("Angle: " + f2.getTheta());
				IJ.log("Phase(psi): " + f2.getPsi());
				IJ.log("Lambda: " + f2.getLambda());
				IJ.log("***");
			} else {
				IJ.log("Filter_Type: Dirac");
				IJ.log("Noise_Level: " + f.getAlpha());
				IJ.log("***");
			}
		}
	}

	/**
	 * this function read a text file and get the plugin's parameters from it
	 * 
	 * @return false if an error is detected in the text file, else return true
	 */
	public boolean readFile() {
		OpenDialog od = new OpenDialog("Choose the file to read", "");

		String path = od.getDirectory() + od.getFileName();
		try {
			// this scanner read the lines
			Scanner scanFile = new Scanner(new File(path));

			// while the file is not end
			while (scanFile.hasNextLine() && this.errorInTxtFile == false) {
				// this scanner read the words
				Scanner scanLine = new Scanner(scanFile.nextLine());
				switch (scanLine.next()) {

				case "Denoising_Algorithm:":
					this.algorithm = scanLine.next().toString();
					// if error
					if (!(this.algorithm.equals("TV"))
							&& !(this.algorithm.equals("H1"))) {

						this.errorInTxtFile = true;
					}
					break;

				case "Iteration_Number:":
					this.nit = Double.parseDouble(scanLine.next());
					break;

				case "Filter_Type:":
					this.filterType = scanLine.next().toString();
					// if error
					if (!(this.filterType.equals("Dirac"))
							&& !(this.filterType.equals("Gabor"))) {
						this.errorInTxtFile = true;
					}
					break;

				case "Noise_Level:":
					this.noiseLevel = Double.parseDouble(scanLine.next());
					if (this.filterType.equals("Dirac")) {
						this.createDirac();
					}
					break;

				case "Sigma_X:":
					this.gammaX = Double.parseDouble(scanLine.next());
					break;

				case "Sigma_Y:":
					this.gammaY = Double.parseDouble(scanLine.next());
					break;

				case "Angle:":
					this.angle = Double.parseDouble(scanLine.next());
					break;

				case "Phase(psi):":
					this.phase_psi = Double.parseDouble(scanLine.next());
					break;

				case "Lambda:":
					this.frequency = Double.parseDouble(scanLine.next());
					if (this.filterType.equals("Gabor")) {
						this.createGabor();
					}
					break;

				case "***":
					break;
				}
				scanLine.close();
			}
			scanFile.close();

		} catch (Exception e) {
			e.printStackTrace();
			IJ.log("Error : the text file is not conform !");
			return false;
		}
		if (this.errorInTxtFile == true) {
			IJ.log("Error : the text file is not conform !");
			return false;
		}
		return true;
	}

	/**
	 * this function ask the user the way to use the plugin, using the GUI or
	 * the Text File (.txt) mode
	 * 
	 * @return false if the window is closed by the user, else return true
	 */
	public boolean askMethod() {
		GenericDialog g = new GenericDialog("Choose the way to use the plugin");

		// signature
		String authors = "Welcome to VSNR plugin !\n \nIn case you use this algorithm, please cite:\nJerome Fehrenbach, Pierre Weiss, Corinne Lorenzo - ITAV\nPlugin Developer : Benjamin Font, Leo Mouly\n \n";
		g.addMessage(authors, new Font(authors,Font.CENTER_BASELINE,13));

		
		String s = "Where are the parameters from ?";
		String[] tabChoice = { "Graphic User Interface (GUI)",
				"Text File (.txt)" };
		String defaultItem = "Choose the input type";
		g.addChoice(s, tabChoice, defaultItem);
		g.pack();
		g.showDialog();

		this.inputMethod = g.getNextChoice();

		if (g.wasCanceled()) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 * thus function ask the user the filter type he want to use
	 * 
	 * @return false if the window is closed by the user, else return true
	 */
	public boolean askFilter() {
		GenericDialog g = new GenericDialog("Filter to use !");

		String[] filterChoice = { "Dirac", "Gabor" };
		g.addChoice("Filter type ?", filterChoice, "Filter...");
		g.pack();
		g.showDialog();

		this.filterType = g.getNextChoice();

		if (g.wasCanceled()) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 * this function ask the user the parameters of the filter
	 * 
	 * @return false if the window is closed by the user, else return true
	 */
	public boolean askParameters() {
		// the filter is dirac
		if (this.filterType.equals("Dirac")) {
			GenericDialog g = new GenericDialog(
					"Setting the Dirac's parameters !");

			g.addNumericField("Noise level :", this.noiseLevel, 0);
			g.pack();
			g.showDialog();

			this.noiseLevel = g.getNextNumber();

			if (g.wasCanceled()) {
				return false;
			} else {
				createDirac();
				return true;
			}
			// the filter is a gabor
		} else {
			GenericDialog g = new GenericDialog(
					"Setting the Gabor's parameters !");

			g.addNumericField("Noise level :", this.noiseLevel, 0);
			g.addNumericField("Sigma X :", this.gammaX, 0);
			g.addNumericField("Sigma Y :", this.gammaY, 0);
			g.addNumericField("Angle :", this.angle, 0);
			g.addNumericField("Phase(psi) :", this.phase_psi, 0);
			g.addNumericField("Lambda :", this.frequency, 0);
			g.pack();
			g.showDialog();

			this.noiseLevel = g.getNextNumber();
			this.gammaX = g.getNextNumber();
			this.gammaY = g.getNextNumber();
			this.angle = g.getNextNumber();
			this.phase_psi = g.getNextNumber();
			this.frequency = g.getNextNumber();

			if (g.wasCanceled()) {
				return false;
			} else {
				createGabor();
				return true;
			}

		}
	}

	/**
	 * this function ask the user if he want to set another filter
	 * 
	 * @return false if the window is closed by the user, else return true
	 */
	public boolean askOtherFilter() {
		GenericDialog g = new GenericDialog("Add another filter ?");
		g.addMessage("Do you want to add an other filter ?");

		g.pack();
		g.showDialog();

		if (g.wasCanceled()) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 * this function ask the user the algorithm to use and the number of
	 * iterations
	 * 
	 * @return false if the window is closed by the user, else return true
	 */
	public boolean askIterationsAndAlgorithm() {
		GenericDialog g = new GenericDialog(
				"Number of iterations and Algorithm ?");

		String[] denoisingChoice = { "TV", "H1" };
		g.addChoice("Choose the denoising algorithm", denoisingChoice, "TV");
		g.addNumericField("Iterations (TV only) :", this.nit, 0);
		IterationListener_2D il = new IterationListener_2D();
		g.addDialogListener(il);
		g.pack();
		g.showDialog();
		
		this.algorithm = g.getNextChoice();
		this.nit = g.getNextNumber();

		if (g.wasCanceled()) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 * this function create a dirac filter using the parameters entered by the
	 * user
	 */
	public void createDirac() {

		Filter_2D dirac = new Filter_2D(FilterType_2D.DIRAC, this.width, this.height, 1);
		dirac.setAlpha(this.noiseLevel);
		// we add the filter to the array
		this.filterList.add(dirac);

	}

	/**
	 * this function create a gabor filter using the parameters entered by the
	 * user
	 */
	public void createGabor() {
		Filter_2D gabor = new GaborFilter_2D(IJ.createImage("Gabor Filter", "32-bit",
				this.image.getWidth(), this.image.getHeight(), 1),
				FilterType_2D.GABOR, this.width, this.height, 1, this.angle,
				this.frequency, this.phase_psi, this.gammaX, this.gammaY);
		gabor.setFlt(new Double2DArray_2D(this.height, this.width));
		gabor.setAlpha(this.noiseLevel);
		// we add the filter to the array
		this.filterList.add(gabor);
	}

	/**
	 * this function denoise the image using the Vsnr_denoiser.denoiseTV_2D()
	 * function
	 * 
	 * @return the denoised image
	 */
	public ImagePlus denoiseTV_2D() {
		// security
		if (filterList.isEmpty()) {

			IJ.log("Unable to process (no filters set)!\nTry to add some filters first \n");
		}
		if (this.image == null) {
			IJ.log("Something bad happened. You probably closed the image.\nPlease reload VSNR! \n");
		}

		// preparation of the parameters of the denoising algorithm
		Double2DArray_2D[] Psis = new Double2DArray_2D[filterList.size()];
		double[] alphas = new double[filterList.size()];

		// we prepare the differents images (Double2DArray)
		Double2DArray_2D img2DToDenoise = new Double2DArray_2D(this.height,
				this.width);
		Double2DArray_2D denoised = new Double2DArray_2D(this.height, this.width);
		Double2DArray_2D noise = new Double2DArray_2D(this.height, this.width);

		// we get the pixels of the image inside a Double2DArray
		ImageProcessor sliceProcessor = this.image.getProcessor();
		for (int i = 0; i < this.height; i++) {
			for (int j = 0; j < this.width; j++) {
				img2DToDenoise.setValue(sliceProcessor.getPixel(j, i), i, j,
						false);
			}
		}

		// this int will help for filling alphas and psis
		int position = 0;

		// for each filter
		for (Filter_2D filter_2D : this.filterList) {
			// if the filter is a gabor, we use the method computeNewGabor()
			if (filter_2D instanceof GaborFilter_2D) {
				GaborFilter_2D gabor = (GaborFilter_2D) filter_2D;
				gabor.computeNewGabor();
				alphas[position] = gabor.getAlpha();
				Psis[position] = gabor.getFlt();
			} else {
				alphas[position] = filter_2D.getAlpha();
				Psis[position] = filter_2D.getFlt();
			}
			position++;
		}

		// we use the denoising algorithm
		denoised = Vsnr_denoiser_2D.denoiseTV_2D(img2DToDenoise, Psis, alphas,
				(int) this.nit, noise, 1, 1);

		// Pierre's code
		Double3DArray_2D denoised3D = new Double3DArray_2D(1, height, width);
		Double3DArray_2D noise3D = new Double3DArray_2D(1, height, width);

		for (int i = 0; i < denoised3D.getRows(); i++) {
			for (int j = 0; j < denoised3D.getColumns(); j++) {
				denoised3D.setValue(denoised.getValue(i, j, false), 0, i, j,
						false);
				noise3D.setValue(noise.getValue(i, j, false), 0, i, j, false);
			}
		}

		ImagePlus resultImg = IJ.createImage("vsnr_result_"
				+ Vsnr_Plugin_2D.this.image.getTitle(), this.type, width,
				height, 1);
		ImageStack resultStack = IJ.createImage(
				"vsnr_result_" + Vsnr_Plugin_2D.this.image.getTitle(),
				this.type, width, height, 1).getImageStack();
		resultStack.deleteLastSlice();

		ImagePlus resultNoise = IJ.createImage("vsnr_noise_"
				+ Vsnr_Plugin_2D.this.image.getTitle(), this.type, width,
				height, 1);
		ImageStack resultStackNoise = IJ.createImage(
				"vsnr_noise_" + Vsnr_Plugin_2D.this.image.getTitle(), this.type,
				width, height, 1).getImageStack();
		resultStackNoise.deleteLastSlice();

		resultImg = IJ.createImage(
				"vsnr_result_" + Vsnr_Plugin_2D.this.image.getTitle(),
				this.type, width, height, 1);
		resultNoise = IJ.createImage(
				"vsnr_noise_" + Vsnr_Plugin_2D.this.image.getTitle(), this.type,
				width, height, 1);

		/*
 		double[] minnmax1 = new double[2];
		double[] minnmax2 = new double[2];
		minnmax1 = denoised3D.minNmax();
		minnmax2 = noise3D.minNmax();
 		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				double val = (denoised3D.getValue(0, i, j, false) - minnmax1[0])
						* 65535.0D / (minnmax1[1] - minnmax1[0]);
				resultImg.getProcessor().putPixel(j, i, (int) Math.floor(val));
				val = (noise3D.getValue(0, i, j, false) - minnmax2[0])
						* 65535.0D / (minnmax2[1] - minnmax2[0]);
				resultNoise.getProcessor()
						.putPixel(j, i, (int) Math.floor(val));
			}
		}

		*/

 		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				double val = denoised3D.getValue(0, i, j, false);
				resultImg.getProcessor().putPixel(j, i, (int) Math.floor(val));
				val = noise3D.getValue(0, i, j, false);
				resultNoise.getProcessor()
						.putPixel(j, i, (int) Math.floor(val));
			}
		}

		
		resultStack.addSlice("Slice #" + 0, resultImg.getProcessor());
		resultStackNoise.addSlice("Slice noise #" + 0,
				resultNoise.getProcessor());

		ImagePlus imresult = new ImagePlus("Denoised Image", resultStack);
		//imresult.setDisplayRange(0.0D, 65535.0D);
		ImagePlus imnoise = new ImagePlus("Noise", resultStackNoise);
		//imnoise.setDisplayRange(0.0D, 65535.0D);
		
		this.noiseImage = imnoise;
		
		this.noiseImage.setDisplayRange(this.image.getDisplayRangeMin(),this.image.getDisplayRangeMax());
		imresult.setDisplayRange(this.image.getDisplayRangeMin(),this.image.getDisplayRangeMax());

		
		return imresult;

	}

	/**
	 * this function denoise the image using the Vsnr_denoiser.denoiseH1_2D()
	 * function
	 * 
	 * @return the denoised image
	 */
	public ImagePlus denoiseH1_2D() {
		// security
		if (filterList.isEmpty()) {

			IJ.log("Unable to process (no filters set)!\nTry to add some filters first \n");
		}
		if (this.image == null) {
			IJ.log("Something bad happened. You probably closed the image.\nPlease reload VSNR! \n");
		}
		// preparation of the parameters of the denoising algorithm
		Double2DArray_2D[] Psis = new Double2DArray_2D[filterList.size()];
		double[] alphas = new double[filterList.size()];

		// we prepare the differents images (Double2DArray)
		Double2DArray_2D img2DToDenoise = new Double2DArray_2D(this.height,
				this.width);
		Double2DArray_2D denoised = new Double2DArray_2D(this.height, this.width);
		Double2DArray_2D noise = new Double2DArray_2D(this.height, this.width);

		// we get the pixels of the image inside a Double2DArray
		ImageProcessor sliceProcessor = this.image.getProcessor();
		for (int i = 0; i < this.height; i++) {
			for (int j = 0; j < this.width; j++) {
				img2DToDenoise.setValue(sliceProcessor.getPixel(j, i), i, j,
						false);
			}
		}

		// this int will help for filling alphas and psis
		int position = 0;

		// for each filter
		for (Filter_2D filter_2D : this.filterList) {
			// if the filter is a gabor, we use the method computeNewGabor()
			if (filter_2D instanceof GaborFilter_2D) {
				GaborFilter_2D gabor = (GaborFilter_2D) filter_2D;
				gabor.computeNewGabor();
				alphas[position] = gabor.getAlpha();
				Psis[position] = gabor.getFlt();
			} else {
				alphas[position] = filter_2D.getAlpha();
				Psis[position] = filter_2D.getFlt();
			}
			position++;
		}

		// we use the denoising algorithm
		denoised = Vsnr_denoiser_2D.denoiseH1_2D(img2DToDenoise, Psis, alphas,
				noise);

		// Pierre's code
		Double3DArray_2D denoised3D = new Double3DArray_2D(1, height, width);
		Double3DArray_2D noise3D = new Double3DArray_2D(1, height, width);

		for (int i = 0; i < denoised3D.getRows(); i++) {
			for (int j = 0; j < denoised3D.getColumns(); j++) {
				denoised3D.setValue(denoised.getValue(i, j, false), 0, i, j,
						false);
				noise3D.setValue(noise.getValue(i, j, false), 0, i, j, false);
			}
		}

		ImagePlus resultImg = IJ.createImage("vsnr_result_"
				+ Vsnr_Plugin_2D.this.image.getTitle(), this.type, width,
				height, 1);
		ImageStack resultStack = IJ.createImage(
				"vsnr_result_" + Vsnr_Plugin_2D.this.image.getTitle(),
				this.type, width, height, 1).getImageStack();
		resultStack.deleteLastSlice();

		ImagePlus resultNoise = IJ.createImage("vsnr_noise_"
				+ Vsnr_Plugin_2D.this.image.getTitle(), this.type, width,
				height, 1);
		ImageStack resultStackNoise = IJ.createImage(
				"vsnr_noise_" + Vsnr_Plugin_2D.this.image.getTitle(), this.type,
				width, height, 1).getImageStack();
		resultStackNoise.deleteLastSlice();


		resultImg = IJ.createImage(
				"vsnr_result_" + Vsnr_Plugin_2D.this.image.getTitle(),
				this.type, width, height, 1);
		resultNoise = IJ.createImage(
				"vsnr_noise_" + Vsnr_Plugin_2D.this.image.getTitle(), this.type,
				width, height, 1);

		/*
		double[] minnmax1 = new double[2];
		
		double[] minnmax2 = new double[2];
		minnmax1 = denoised3D.minNmax();
		minnmax2 = noise3D.minNmax();
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				double val = (denoised3D.getValue(0, i, j, false) - minnmax1[0])
						* 65535.0D / (minnmax1[1] - minnmax1[0]);
				resultImg.getProcessor().putPixel(j, i, (int) Math.floor(val));
				val = (noise3D.getValue(0, i, j, false) - minnmax2[0])
						* 65535.0D / (minnmax2[1] - minnmax2[0]);
				resultNoise.getProcessor()
						.putPixel(j, i, (int) Math.floor(val));
			}
		}
		*/
		
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				double val = denoised3D.getValue(0, i, j, false);
				resultImg.getProcessor().putPixel(j, i, (int) Math.floor(val));
				val = noise3D.getValue(0, i, j, false);
				resultNoise.getProcessor()
						.putPixel(j, i, (int) Math.floor(val));
			}
		}
		
		resultStack.addSlice("Slice #" + 0, resultImg.getProcessor());
		resultStackNoise.addSlice("Slice noise #" + 0,
				resultNoise.getProcessor());

		ImagePlus imresult = new ImagePlus("Denoised Image", resultStack);
		//imresult.setDisplayRange(0.0D, 65535.0D);
		ImagePlus imnoise = new ImagePlus("Noise", resultStackNoise);
		//imnoise.setDisplayRange(0.0D, 65535.0D);
		this.noiseImage = imnoise;
		
		this.noiseImage.setDisplayRange(this.image.getDisplayRangeMin(),this.image.getDisplayRangeMax());
		imresult.setDisplayRange(this.image.getDisplayRangeMin(),this.image.getDisplayRangeMax());
		
		return imresult;
	}

}
