package views;

import ij.ImageJ;

public class Test_2D {
	public static void main(String[] args) {
		Class<?> c = Vsnr_Plugin_2D.class;
		String url = c.getResource(
				"/" + c.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length()
				- c.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		
		// Launch Fiji
		new ImageJ();
	}
}