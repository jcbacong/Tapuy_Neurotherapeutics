/*
 * Macro template to process multiple images in a folder
 */

#@File(label = "Input directory", style = "directory") input
#@File (label = "Output directory", style = "directory") output
#@String (label = "File suffix", value = ".png") suffix

// See also Process_Folder.py for a version of this code
// in the Python scripting language.

setBatchMode(true); // Prevent opening new windows
processFolder(input);
print("Done batch processing...");

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	setBatchMode(true);
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i], suffix);
	}
}

function processFile(input, output, file, suffix) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	setBatchMode(true);
	figfile = input + File.separator + file;
	scale = "/Users/jcbacong/Documents/BCH 205/BCH205-Tapuy-Project/Retardation Assay/scale.png";
	
	// Remove the scale and bg noise
	open(figfile);
	open(scale);
	selectWindow(file);
	imageCalculator("Subtract create", file, "scale.png");
	selectWindow("Result of "+file);
	
	// Image processing
	// Median filtering of GFP
	run("Split Channels");
	selectWindow("Result of "+ file+" (green)");
	run("Despeckle");
	run("Duplicate...", "title=median-gfp");
	run("Median...", "radius=2");
	imageCalculator("Subtract create", "Result of "+ file+" (green)","median-gfp");
	selectWindow("Result of Result of "+file+" (green)");
	
	// Background subtraction & Thresholding
	run("Despeckle");
	run("Subtract Background...", "rolling=2.0");


	run("Auto Threshold", "method=RenyiEntropy white");


	// Generating ROIs and measurement
	run("Set Measurements...", "area mean standard min integrated redirect=[Result of "+file+" (green)] decimal=3");
	run("Analyze Particles...", "size=2-150 circularity=0.01-1.00 show=[Overlay Masks] display exclude add");
	
	// Save results from table
	saveAs("Results", output+ File.separator + replace(file,suffix,"") + " Results.csv");
	
	run("Flatten");
	// roiManager("Show All without labels");
	// Save neuron-highlighted images
	saveAs("Tiff", output+ File.separator + replace(file,suffix,"") + " Output");
	run("Close");
	
	print("Processing: " + input + File.separator + file);
	print("Saving to: " + output);
	
	// Close windows

	selectWindow("Result of "+ file+" (red)");
	run("Close");
	selectWindow("Result of "+ file+" (blue)");
	run("Close");
	selectWindow("Result of "+ file+" (green)");
	run("Close");
	selectWindow("median-gfp");
	run("Close");
	selectWindow("Result of Result of "+file+" (green)");
	run("Close");
	
}
