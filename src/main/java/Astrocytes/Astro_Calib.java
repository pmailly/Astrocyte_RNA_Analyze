package Astrocytes;



import static Tools.Astro_Dots_Tools.bg;
import static Tools.Astro_Dots_Tools.dialog;
import static Tools.Astro_Dots_Tools.find_background;
import static Tools.Astro_Dots_Tools.find_dots;
import static Tools.Astro_Dots_Tools.flush_close;
import static Tools.Astro_Dots_Tools.labelsObject;
import static Tools.Astro_Dots_Tools.localThickness3D;
import static Tools.Astro_Dots_Tools.objectsSizeFilter;
import static Tools.Astro_Dots_Tools.stdBg;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
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
import mcib3d.image3d.ImageHandler;
import mcib3d.utils.ArrayUtil;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */
        
        
public class Astro_Calib implements PlugIn {
    
    private boolean canceled;
    private String imageDir;
    private String outDirResults;
    private Calibration cal = new Calibration();
    private BufferedWriter outPutResults;
    private final double minDotsSize = 0.03;
    private final double maxDotsSize = 20;
    private final String[] ch = {"0", "1", "2", "3"}; 
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
            imageDir = IJ.getDirectory("Choose Directory Containing ICS Files...");
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
                if (imageFile[i].endsWith(".ics")) {
                    String imageName = inDir+ File.separator+imageFile[i];
                    String rootName = imageFile[i].replace(".ics", "");
                    // Find ROI calibration file
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
                            // return the index for channels 0 DAPI, 1 Astro, 2 Dots and ask for calibration if needed 
                            chIndex = dialog(ch, showCal, cal);
                            if (ch == null)
                                return;
                            cal.setUnit("microns"); 
                            System.out.println("x cal = " +cal.pixelWidth+", z cal = " + cal.pixelDepth);
                        }
                        // Write headers for results file
                        FileWriter fileResults = new FileWriter(outDirResults + rootName+"_calibration_results.xls", false);
                        outPutResults = new BufferedWriter(fileResults);
                        outPutResults.write("Roi Name\tBackground\tStd background\tMin Intensity\tMax Intensity\tMean Vol dot\tSdt Vol\n");
                        outPutResults.flush();
                        // find rois
                        RoiManager rm = new RoiManager(false);
                        rm.runCommand("Open", roi_file);
                        // for each roi open image and crop
                        
                        // Write headers for results file
                        FileWriter fileRoiResults = new FileWriter(outDirResults + rootName + "_RoiResults.xls", false);
                        BufferedWriter outRoiResults = new BufferedWriter(fileRoiResults);
                        outRoiResults.write("Roi Name\t#Dot\tDot mean intensity in Astro\tDot Z\tDot Vol\tAstro diameter\n");
                        fileRoiResults.flush();
                        int roiIndex = 0;
                        ImporterOptions options = new ImporterOptions();
                        options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                        options.setId(imageName);
                        options.setSplitChannels(true);
                        
                        /** 
                         * open channels
                         */
                        // astro channel    
                        ImagePlus imgAstro = BF.openImagePlus(options)[chIndex.get(1)];
                        IJ.showStatus("Opening Astrocyte channel");
                       
                        // dots channel
                        ImagePlus imgDots = BF.openImagePlus(options)[chIndex.get(2)];
                        IJ.showStatus("Opening Dots channel");

                        for (int r = 0; r < rm.getCount(); r++) {
                            roiIndex++;
                            rm.select(imgAstro, r);
                            imgAstro.updateAndDraw();
                            Roi roi = imgAstro.getRoi();
                            
                            
                            // Find in roi name the desired top and bottom stack 
                            // roi name should be roi_number-ztop-zbottom
                            String[] regExp = rm.getName(r).split("-");
                            if (Integer.parseInt(regExp[1]) < 1)
                                regExp[1] = "1";
                            if (Integer.parseInt(regExp[2]) > imgAstro.getNSlices())
                                regExp[2] = Integer.toString(imgAstro.getNSlices());
                            int zStart = Integer.parseInt(regExp[1]);
                            int zStop = Integer.parseInt(regExp[2]);
                            
                            // make substack and crop with roi if exist                           
                            ImagePlus imgAstroZCrop = new Duplicator().run(imgAstro, zStart, zStop);
                            imgAstroZCrop.setTitle(rootName+"_Astro");
                            imgAstroZCrop.deleteRoi();
 
                            // Compute distance map
                            ImagePlus imgAstroZCropMap = localThickness3D(imgAstroZCrop);
                            IJ.run(imgAstroZCropMap, "Fire", "");
                            imgAstroZCropMap.setCalibration(cal);
                            FileSaver imgAstroMapFile = new FileSaver(imgAstroZCropMap);
                            imgAstroMapFile.saveAsTiff(outDirResults + rootName + "_Astro-"+roiIndex+"-map.tif"); 
                            
                            // open dots channel
                            rm.select(imgDots, r);
                            imgDots.updateAndDraw();
                            ImagePlus imgDotsZCrop = new Duplicator().run(imgDots, zStart, zStop);

                            // Find dots population
                            Objects3DPopulation dotsPop = find_dots(imgDotsZCrop, imgDotsZCrop.getRoi());
                            System.out.println("Dots number = "+dotsPop.getNbObjects());
                            objectsSizeFilter(minDotsSize, maxDotsSize, dotsPop, imgDotsZCrop, false);
                            System.out.println("After size filter dots number = "+dotsPop.getNbObjects());
                            
                            // save dots image
                            ImageHandler imgObjDots = ImageHandler.wrap(imgDotsZCrop).createSameDimensions();
                            for (int n = 0; n < dotsPop.getNbObjects(); n++) {
                                Object3D obj = dotsPop.getObject(n);
                                obj.draw(imgObjDots, 255);
                                labelsObject(obj, imgObjDots.getImagePlus(), n, 255, 10);
                            }
                            
                            
                            // compute min and dots Intensity
                            ArrayUtil statsIntensity = new ArrayUtil(dotsPop.getNbObjects());
                            ArrayUtil statsVolume = new ArrayUtil(dotsPop.getNbObjects());
                            ImageHandler imgH = ImageHandler.wrap(imgAstroZCrop);
                            double meanInt = 0;
                            double volObj = 0; 
                            double astroDiameter = 0;
                            for (int n = 0; n < dotsPop.getNbObjects(); n++) {
                                Object3D obj = dotsPop.getObject(n);
                                volObj = obj.getVolumeUnit();
                                astroDiameter = obj.getPixMaxValue(ImageHandler.wrap(imgAstroZCropMap)) * cal.pixelWidth;
                                meanInt = obj.getPixMeanValue(imgH);
                                outRoiResults.write(roi.getName()+"\t"+n+"\t"+meanInt+"\t"+obj.getCenterZ()+"\t"+volObj+"\t"+astroDiameter+"\n");
                                outRoiResults.flush();
                                statsIntensity.addValue(n, meanInt); 
                                statsVolume.addValue(n, volObj);
                            }
                            double minInt = statsIntensity.getMinimum();
                            double maxInt = statsIntensity.getMaximum();
                            find_background(imgAstroZCrop);
                            double meanVol = statsVolume.getMean();
                            double stdVol = statsVolume.getStdDev();
                            
                            outPutResults.write(roi.getName()+"\t"+bg+"\t"+stdBg+"\t"+minInt+"\t"+maxInt+"\t"+meanVol+"\t"+stdVol+"\n");
                            
                            // save image for objects population
                            ImagePlus[] imgColors = {imgObjDots.getImagePlus(), null, null, imgAstroZCrop};
                            ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
                            imgObjects.setCalibration(imgAstroZCrop.getCalibration());
                            IJ.run(imgObjects, "Enhance Contrast", "saturated=0.35");
                            FileSaver imgObjectsFile = new FileSaver(imgObjects);
                            imgObjectsFile.saveAsTiff(outDirResults + rootName + "_Astro-"+roiIndex+"-Dots.tif"); 
                            flush_close(imgObjects);
                            imgAstroZCrop.close();
                            imgDotsZCrop.close();
                        }
                        flush_close(imgAstro);
                        flush_close(imgDots);
                        outRoiResults.close();
                    }
                } 
            }
            outPutResults.close();
            IJ.showStatus("Calibration done....");
        } catch (DependencyException | ServiceException | IOException | FormatException ex) {
            Logger.getLogger(Astro_Calib.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}


