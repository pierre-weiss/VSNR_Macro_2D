package views;

import java.awt.AWTEvent;
import java.awt.Choice;
import java.awt.TextField;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;

public class IterationListener_2D implements DialogListener {

	@Override
	public boolean dialogItemChanged(GenericDialog g, AWTEvent e) {
		if (((Choice) g.getChoices().get(0)).getSelectedItem().equals("TV")) {
			((TextField) g.getNumericFields().get(0)).setEnabled(true);
			return true;
		} else {
			((TextField) g.getNumericFields().get(0)).setEnabled(false);
			return true;
		}
	}
}