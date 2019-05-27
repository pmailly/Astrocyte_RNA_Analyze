package Astrocytes;


import static Tools.Astro_Dots_Tools.*;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */



public class Astro_Dots implements PlugIn {
    private static String imageDir;
    private final boolean canceled = false;
    public static String outDirResults = "";
    private Calibration cal = new Calibration();
    private final double stepZ = 0.3;
    private final double minNucSize = 50;
    private final double maxNucSize = 5000;
    private final double minDotsSize = 0.03;
    private final double maxDotsSize = 2000;
    public static double meanDotsSize = 0.6;
    public static double stdMeanDotsSize = 0.6;
    private BufferedWriter outPutResults;
    
    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir = IJ.getDirectory("Choose Directory Containing ND Files...");
            if (imageDir == null) {
                return;
            }
            File inDir = new File(imageDir);
            String[] imageFile = inDir.list();
            if (imageFile == null) {
                return;
            }
            // create output folder
            outDirResults = inDir + File.separator+ "Out"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            // Write headers for results file
            FileWriter fileResults = new FileWriter(outDirResults + "Astro_results.xls", false);
            outPutResults = new BufferedWriter(fileResults);
            outPutResults.write("ImageName\tRoi Name\tMean background\tStd background\tAstrocyte Volume\tDensity dots in Astro"
                    + "\tPercentage of dots not in astro\tPercentage of dots in soma\tPercentage of dots in fine processes\tPercentage of dots in large processes"
                    + "\tDots mean intensity in Astro\tSD intensity in astro\tMean astro diameter(0 exluded)\tStd astro diameter(0 excluded)"
                    + "\tMed astro diameter\tMean nucleus distance\tStd nucleus distance\n");
            outPutResults.flush();
            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            Arrays.sort(imageFile);
            int imageNum = 0;
            ArrayList<Integer> chIndex = new ArrayList();
            for (int i = 0; i < imageFile.length; i++) {
                // Find nd files
                if (imageFile[i].endsWith(".nd")) {
                    String imageName = inDir+ File.separator+imageFile[i];
                    String rootName = imageFile[i].replace(".nd", "");
                    // Find ROI file
                    String roi_file = imageDir+rootName+".zip";
                    if (!new File(roi_file).exists()) {
                        IJ.showStatus("No ROI file found !");
                        return;
                       }
                    else {
                        reader.setId(imageName);
                        int series = 0;
                        reader.setSeries(0);
                        imageNum++;
                        boolean showCal = false;
                        // Check calibration
                        if (imageNum == 1) {
                            cal.pixelWidth = meta.getPixelsPhysicalSizeX(series).value().doubleValue();
                            cal.pixelHeight = cal.pixelWidth;
                            // problem to read calibration with nd files
                            if ((meta.getPixelsPhysicalSizeZ(series) == null) || (cal.pixelWidth == 1))
                                showCal = true;
                            else
                                cal.pixelDepth = meta.getPixelsPhysicalSizeZ(series).value().doubleValue();
                            String[] seriesName = meta.getImageName(0).split("/");
                            // return the index for channels 0 DAPI, 1 Astro, 2 Dots and ask for calibration if needed 
                            chIndex = dialog(seriesName, showCal, cal, true);
                            if (chIndex == null)
                                return;
                            cal.setUnit("microns");
                            System.out.println("x cal = " +cal.pixelWidth+", z cal = " + cal.pixelDepth);
                        }
                        // find rois
                        RoiManager rm = new RoiManager(false);
                        rm.runCommand("Open", roi_file);
                        int index = 0;
                        ImporterOptions options = new ImporterOptions();
                        options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                        options.setId(imageDir + File.separator + rootName + "_"+ rootName + "_0_cmle.ics");
                        options.setSplitChannels(true);
                        /**
                         * Open channels
                         */
                        // Nucleus channel
                        ImagePlus imgNuc = new ImagePlus();  
                        if (useDeconvImg) {
                            IJ.showStatus("Opening Nucleus channel");
                            imgNuc = BF.openImagePlus(options)[chIndex.get(0)];
                        }
                        else
                            imgNuc = spinningReadChannel(reader, chIndex.get(0), imageName, cal);
                        
                        // Astrocyte channel
                        ImagePlus imgAstro = new ImagePlus();  
                        if (useDeconvImg) {
                            IJ.showStatus("Opening Astrocyte channel");
                            imgAstro = BF.openImagePlus(options)[chIndex.get(1)];
                        }
                        else
                            imgAstro = spinningReadChannel(reader, chIndex.get(1), imageName, cal);
                        
                        // Dots channel
                        ImagePlus imgDots = new ImagePlus();  
                        if (useDeconvImg) {
                            IJ.showStatus("Opening Dots channel");
                            imgDots = BF.openImagePlus(options)[chIndex.get(2)];
                        }
                        else
                            imgDots = spinningReadChannel(reader, chIndex.get(2), imageName, cal);
                        
                        // for each roi open image and crop
                        for (int r = 0; r < rm.getCount(); r++) {
                            index++;                            
                            // Find in roi name the desired top and bottom stack 
                            // roi name should be roi_number-ztop-zbottom
                            String[] regExp = rm.getName(r).split("-");
                            if (Integer.parseInt(regExp[1]) < 1)
                                regExp[1] = "1";
                            if (Integer.parseInt(regExp[2]) > imgNuc.getNSlices())
                                regExp[2] = Integer.toString(imgNuc.getNSlices());
                            int zStart = Integer.parseInt(regExp[1]);
                            int zStop = Integer.parseInt(regExp[2]);
                            rm.select(imgNuc,r);
                            imgNuc.updateAndDraw();
                           
                            // make substack
                            ImagePlus imgNucZCrop = new Duplicator().run(imgNuc, zStart, zStop);
                            rm.select(imgNucZCrop,r);
                            Roi roi = imgNucZCrop.getRoi();
                            imgNucZCrop.deleteRoi();
                            median_filter(imgNuc, 2);
                            IJ.run(imgNucZCrop, "Difference of Gaussians", "  sigma1=30 sigma2=20 stack");
                            threshold(imgNucZCrop, "Li", true);
                            imgNucZCrop.setRoi(roi);
                            imgNucZCrop.updateAndDraw();
                            IJ.run("Colors...", "foreground=black background=white selection=yellow");
                            IJ.run(imgNucZCrop, "Clear Outside", "stack");
                            imgNucZCrop.deleteRoi();
                            IJ.run(imgNucZCrop, "Options...", "iterations=10 count=1 do=Open stack");
                            IJ.run(imgNucZCrop, "Fill Holes", "stack");
                            
                            // WaterShed slipt
                            ImagePlus imgNucSplit = watershedSplit(imgNucZCrop, 6);
                            imgNucSplit.setCalibration(cal);
                            IJ.run(imgNucSplit, "Fill Holes", "stack");
                            // find nucleus population
                            Objects3DPopulation nucPop = getPopFromImage(imgNucSplit, cal);
                            System.out.println("Roi = "+index+" of " + rm.getCount());
                            System.out.println("Nucleus number = "+nucPop.getNbObjects());
                            objectsSizeFilter(minNucSize, maxNucSize, nucPop, imgNucSplit, true);
                            System.out.println("After size filter Nucleus number = "+nucPop.getNbObjects());
                            flush_close(imgNucZCrop);
                            flush_close(imgNucSplit);
                            if (nucPop.getNbObjects() != 0) {  
                                
                                // astro channel  
                                rm.select(imgAstro, r);
                                imgAstro.updateAndDraw();
                                // make substack
                                ImagePlus imgAstroZCrop = new Duplicator().run(imgAstro, zStart, zStop);  
                                imgAstroZCrop.setTitle(rootName+"_Astro");
                                
                                Object3D nucAstro = null;
                                // if more than one nucleus find nucleus with intensity in astrocyte channel
                                // and ask to choose
                                if (nucPop.getNbObjects() > 1) {
                                    ImagePlus img = imgAstroZCrop.duplicate();
                                    nucAstro = nucleusSelect(img, nucPop);
                                    flush_close(img);
                                }                                    
                                else
                                    nucAstro = nucPop.getObject(0);
                                
                                // compute distance map image
                                
                                ImagePlus imgAstroZCropMap = localThickness3D(imgAstroZCrop);
                                imgAstroZCropMap.setCalibration(cal);
                                
                                // dots channel
                                rm.select(imgDots, r);
                                imgDots.updateAndDraw();
                                ImagePlus imgDotsZCrop = new Duplicator().run(imgDots, zStart, zStop); 
                                
                                // Find dots population
                                Objects3DPopulation dotsPop = find_dots(imgDotsZCrop, roi, useDeconvImg);
                                System.out.println("Dots number = "+dotsPop.getNbObjects());
                                objectsSizeFilter(minDotsSize, maxDotsSize, dotsPop, imgDotsZCrop, false);
                                System.out.println("After size filter dots number = "+dotsPop.getNbObjects());
                                
                                // calculate parameters
                                classify_dots(nucAstro, dotsPop, imgAstroZCrop, imgAstroZCropMap);                              

                                // draw objects
                                tagsObjects(nucAstro, dotsPop, imgAstroZCrop, outDirResults, rootName, r);                               

                                // write a global parameters table by image
                                compute_Image_parameters(roi, imgAstroZCrop, imgAstroZCropMap, nucAstro, dotsPop, outPutResults, rootName);
                                flush_close(imgAstroZCrop);
                                flush_close(imgAstroZCropMap);
                            }
                            else 
                                System.out.println("No nucleus found !");
                        }
                        
                        flush_close(imgNuc);
                        flush_close(imgAstro);
                        flush_close(imgDots);
                    }
                }
            }
            outPutResults.close();
        }
        catch (IOException | FormatException | DependencyException | ServiceException ex) {
            Logger.getLogger(Astro_Dots.class.getName()).log(Level.SEVERE, null, ex);
        }
        IJ.showStatus("Process done ...");
    }
}
