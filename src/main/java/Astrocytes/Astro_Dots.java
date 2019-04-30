package Astrocytes;


import static Tools.Astro_Dots_Tools.*;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.plugin.SubstackMaker;
import ij.plugin.frame.RoiManager;
import ij.process.AutoThresholder;
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
    private final double minNucSize = 100;
    private final double maxNucSize = 5000;
    private final double minDotsSize = 0.2;
    private final double maxDotsSize = 20;
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
                            cal.setUnit("microns");
                            System.out.println("x cal = " +cal.pixelWidth+", z cal = " + cal.pixelDepth);
                        }
                        // find rois
                        RoiManager rm = new RoiManager(false);
                        rm.runCommand("Open", roi_file);
                        // for each roi open image and crop
                        for (int r = 0; r < rm.getCount(); r++) {
                            // Write headers for results file
                            FileWriter fileResults = new FileWriter(outDirResults + rootName + "_Astro_" + r +"_results.xls", false);
                            outPutResults = new BufferedWriter(fileResults);
                            outPutResults.write("Roi Name\tAstrocyte Volume\tAstrocyte circularity\tNucleus Vol\tSoma Vol\t#Dot\tDot Vol\tDot Norm intDensity in Astro"
                                    + "\tDistance to nucleus border\tDots type\tBg Int\n");
                            outPutResults.flush();
                            
                            // open nucleus channel
                            ImagePlus imgNuc = spinningReadChannel(reader, chIndex.get(0), imageName, cal);
                            rm.select(imgNuc,r);
                            imgNuc.updateAndDraw();
                            IJ.run(imgNuc, "Crop","");
                            Roi roi = imgNuc.getRoi();
                            imgNuc.deleteRoi();
                            // Find in roi name the desired top and bottom stack 
                            // roi name should be roi_number-ztop-zbottom
                            String[] regExp = rm.getName(r).split("-");
                            if (Integer.parseInt(regExp[1]) < 1)
                                regExp[1] = "1";
                            if (Integer.parseInt(regExp[2]) > imgNuc.getNSlices())
                                regExp[2] = Integer.toString(imgNuc.getNSlices());
                            String range = regExp[1] + "-" + regExp[2];
                            // make substack
                            ImagePlus imgNucCrop = new SubstackMaker().makeSubstack(imgNuc, range);
                            IJ.run(imgNucCrop, "Difference of Gaussians", "  sigma1=30 sigma2=15 enhance stack");
                            threshold(imgNucCrop, AutoThresholder.Method.Li, true);
                            if (roi != null) {
                                imgNucCrop.setRoi(roi);
                                imgNucCrop.updateAndDraw();IJ.run("Colors...", "foreground=black background=white selection=yellow");
                                IJ.run(imgNucCrop, "Clear Outside", "stack");
                                imgNucCrop.deleteRoi();
                            }
                            IJ.run(imgNucCrop, "Options...", "iterations=10 count=1 do=Open stack");
                            IJ.run(imgNucCrop, "Fill Holes", "stack");
                            // WaterShed slipt
                            ImagePlus imgNucSplit = watershedSplit(imgNucCrop, 2);
                            imgNucSplit.setCalibration(cal);
                            IJ.run(imgNucSplit, "Fill Holes", "stack");
                            // find nucleus population
                            Objects3DPopulation nucPop = getPopFromImage(imgNucSplit, cal);
                            System.out.println("Roi = "+r+" of " + rm.getCount());
                            System.out.println("Nucleus number = "+nucPop.getNbObjects());
                            objectsSizeFilter(minNucSize, maxNucSize, nucPop, imgNucSplit, true);
                            System.out.println("After size filter Nucleus number = "+nucPop.getNbObjects());
                            imgNuc.close();
                            flush_close(imgNucCrop);
                            flush_close(imgNucSplit);
                            if (nucPop.getNbObjects() != 0) {
                                
                                // Open astro channel
                                ImagePlus imgAstro = spinningReadChannel(reader, chIndex.get(1), imageName, cal);    
                                rm.select(imgAstro, r);
                                imgAstro.updateAndDraw();
                                IJ.run(imgAstro, "Crop","");
                                imgAstro.deleteRoi();
                                // bleach correction
                                if (doBleachCorr) 
                                    doBleachCorrection(imgAstro);
                                // make substack
                                ImagePlus imgAstroZCrop = new SubstackMaker().makeSubstack(imgAstro, range);
                                imgAstro.close();
                                
                                Object3D nucAstro = null;
                                // if more than one nucleus find nucleus with intensity in astrocyte channel
                                // and ask to choose
                                if (nucPop.getNbObjects() > 1) 
                                    nucAstro = nucleusSelect(imgAstroZCrop, nucPop);
                                else
                                    nucAstro = nucPop.getObject(0);
                                // Open dots channel
                                ImagePlus imgDots = spinningReadChannel(reader, chIndex.get(2), imageName, cal);
                                rm.select(imgDots, r);
                                imgDots.updateAndDraw();
                                IJ.run(imgDots, "Crop","");
                                ImagePlus imgDotsCrop = new SubstackMaker().makeSubstack(imgDots, range);
                                imgDots.close();
                                

                                // Find dots population
                                Objects3DPopulation dotsPop = find_dots(imgDotsCrop, roi);
                                System.out.println("Dots number = "+dotsPop.getNbObjects());
                                objectsSizeFilter(minDotsSize, maxDotsSize, dotsPop, imgDotsCrop, false);
                                System.out.println("After size filter dots number = "+dotsPop.getNbObjects());
                                
                                // calculate parameters
                                compute_parameters(roi, nucAstro, dotsPop, imgDotsCrop, imgAstroZCrop, outPutResults);                              
                               
                                // draw objects
                                tagsObjects(nucAstro, dotsPop, imgAstroZCrop, outDirResults, rootName, r);                               
                                flush_close(imgAstroZCrop);
                                flush_close(imgDotsCrop);
                                outPutResults.close();
                            }
                            else {
                                System.out.println("No nucleus found !");
                                outPutResults.close();
                                // delete empty file
                                new File(outDirResults + rootName + "_Astro_" + r +"_results.xls").delete();
                            }    
                        }
                        // write a global parameters table by image
                        compute_Image_parameters(outDirResults, rootName);
                    }
                }
            }
        }
        catch (IOException | FormatException | DependencyException | ServiceException ex) {
            Logger.getLogger(Astro_Dots.class.getName()).log(Level.SEVERE, null, ex);
        }
        IJ.showStatus("Process done ...");
    }
}
