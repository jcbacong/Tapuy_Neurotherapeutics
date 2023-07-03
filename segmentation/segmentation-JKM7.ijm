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
//	blank = "/Users//blank 06-24-2022.tif";
	
	// Remove the scale and bg noise
	open(figfile);
	
	// Image processing
	run("Split Channels");
	selectWindow(file+" (red)");
	run("Duplicate...", "title=source");
	
	selectWindow(file+" (red)");
	run("Sharpen");
	
	// Background subtraction & Thresholding
//	run("Subtract Background...", "rolling=10 sliding disable");
	run("Auto Threshold", "method=Otsu white");
//	run("Auto Threshold", "method=Minimum white");
//	run("Auto Threshold", "method=Shanbhag white");
//	run("Auto Threshold", "method=IsoData white");

	// Generating ROIs and measurement
	run("Set Measurements...", "area mean standard min integrated redirect=source decimal=3");
	run("Analyze Particles...", "size=10-Infinity circularity=0.01-Infinity show=[Overlay Masks] display exclude add");
	// Save results from table
	saveAs("Results", output+ File.separator + replace(file,suffix,"") + " Results.csv");
	
	run("Flatten");
	// roiManager("Show All without labels");
	// Save neuron-highlighted images
	saveAs("Tiff", output+ File.separator + replace(file,suffix,"") + " Output");
	run("Close");
	
	print("Processing: " + input + File.separator + file);
	print("Saving to: " + output);


	selectWindow(file+" (green)");
	close();
	selectWindow(file+" (red)");
	close();
	selectWindow(file+" (blue)");
	close();
	selectWindow("source");
	close();
}
