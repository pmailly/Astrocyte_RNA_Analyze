package Astrocytes;



import static Tools.Astro_Dots_Tools.dialog;
import static Tools.Astro_Dots_Tools.doBleachCorr;
import static Tools.Astro_Dots_Tools.doBleachCorrection;
import static Tools.Astro_Dots_Tools.find_background;
import static Tools.Astro_Dots_Tools.find_dots;
import static Tools.Astro_Dots_Tools.objectsSizeFilter;
import static Tools.Astro_Dots_Tools.spinningReadChannel;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.plugin.SubstackMaker;
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
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

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
    private final double minDotsSize = 0.2;
    private final double maxDotsSize = 20;

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
                    // Find ROI calibration file
                    String roi_file = imageDir+rootName+"_cal.zip";
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
                            chIndex = dialog(seriesName, showCal, cal, false);
                            cal.setUnit("microns");
                            System.out.println("x cal = " +cal.pixelWidth+", z cal = " + cal.pixelDepth);
                        }
                        // Write headers for results file
                        FileWriter fileResults = new FileWriter(outDirResults + "calibration_results.xls", false);
                        outPutResults = new BufferedWriter(fileResults);
                        outPutResults.write("Roi Name\tBackground\tMin Intensity\tMax Intensity\tMin ratio\tMax ratio\n");
                        outPutResults.flush();
                        // find rois
                        RoiManager rm = new RoiManager(false);
                        rm.runCommand("Open", roi_file);
                        // for each roi open image and crop
                        
                        
                        for (int r = 0; r < rm.getCount(); r++) {
                            // open astrocyte channel
                            ImagePlus imgAstro = spinningReadChannel(reader, chIndex.get(1), imageName, cal);
                            rm.select(imgAstro, r);
                            imgAstro.updateAndDraw();
                            Roi roi = imgAstro.getRoi();
                            IJ.run(imgAstro, "Crop","");
                            imgAstro.deleteRoi();
                            // Find in roi name the desired top and bottom stack 
                            // roi name should be roi_number-ztop-zbottom
                            String[] regExp = rm.getName(r).split("-");
                            if (Integer.parseInt(regExp[1]) < 1)
                                regExp[1] = "1";
                            if (Integer.parseInt(regExp[2]) > imgAstro.getNSlices())
                                regExp[2] = Integer.toString(imgAstro.getNSlices());
                            String range = regExp[1] + "-" + regExp[2];
                            
                            // bleach correction
                            if (doBleachCorr)
                                doBleachCorrection(imgAstro);
                            // make substack                            
                            ImagePlus imgAstroZCrop = new SubstackMaker().makeSubstack(imgAstro, range);
                            imgAstro.close();
                            
                            // open dots channel
                            ImagePlus imgDots = spinningReadChannel(reader, chIndex.get(2), imageName, cal);
                            rm.select(imgDots, r);
                            imgDots.updateAndDraw();
                            IJ.run(imgDots, "Crop","");
                            ImagePlus imgDotsCrop = new SubstackMaker().makeSubstack(imgDots, range);
                            imgDots.close();

                            // Find dots population
                            Objects3DPopulation dotsPop = find_dots(imgDotsCrop, imgDots.getRoi());
                            System.out.println("Dots number = "+dotsPop.getNbObjects());
                            objectsSizeFilter(minDotsSize, maxDotsSize, dotsPop, imgDotsCrop, false);
                            System.out.println("After size filter dots number = "+dotsPop.getNbObjects());
                            
                            // compute min and dots Intensity
                            DescriptiveStatistics statsIntensity = new DescriptiveStatistics();
                            ImageHandler imgH = ImageHandler.wrap(imgAstroZCrop);
                            for (int n = 0; n < dotsPop.getNbObjects(); n++) {
                                Object3D obj = dotsPop.getObject(n);
                                statsIntensity.addValue(obj.getIntegratedDensity(imgH)/obj.getVolumeUnit());    
                            }
                            double minInt = statsIntensity.getMin();
                            double maxInt = statsIntensity.getMax();
                            double bg = find_background(imgAstroZCrop);
                            double minRatio = minInt/bg;
                            double maxRatio = maxInt/bg;
                            outPutResults.write(roi.getName()+"\t"+bg+"\t"+minInt+"\t"+maxInt+"\t"+minRatio+"\t"+maxRatio+"\n");
                            imgAstroZCrop.close();
                            imgDotsCrop.close();
                        }
                    }
                }
                outPutResults.close();
                if(doBleachCorr)
                    WindowManager.getWindow("Log").dispose();
                IJ.showStatus("Calibration done....");
            }
        } catch (DependencyException | ServiceException | IOException | FormatException ex) {
            Logger.getLogger(Astro_Calib.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}


