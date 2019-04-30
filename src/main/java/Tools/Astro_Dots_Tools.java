package Tools;

import emblcmci.BleachCorrection_SimpleRatio;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.GaussianBlur3D;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.RankFilters;
import ij.process.AutoThresholder;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.awt.Checkbox;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.formats.FormatException;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import mcib3d.geom.Object3DVoxels;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.image3d.processing.FastFilters3D;
import mcib3d.image3d.regionGrowing.Watershed3D;
import mcib3d.utils.ThreadUtil;
import org.apache.commons.io.filefilter.RegexFileFilter;
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


public class Astro_Dots_Tools {
    
//Intensity threshold for dots in GFAP channel   
public static boolean doBleachCorr = true; 
public static double minThreshold = 0;
public static double maxThreshold = 65000;

    /**
     * Read channels
     * @param r
     * @param channel
     * @param name
     * @param cal
     * @return image
     */
    public static ImagePlus spinningReadChannel(ImageProcessorReader r, int channel, String name, Calibration cal) {
        ImageStack stack = new ImageStack(r.getSizeX(), r.getSizeY());
        int start = r.getSizeZ() * channel;
        int stop = r.getSizeZ() - 1 + r.getSizeZ() * channel;
        IJ.showStatus("reading channel "+channel+" ...");
        for (int n = start; n <= stop; n++) {
            ImageProcessor ip = null;
            try {
                try {
                    ip = r.openProcessors(n)[0];
                } catch (FormatException ex) {
                    Logger.getLogger(Astro_Dots_Tools.class.getName()).log(Level.SEVERE, null, ex);
                }
            } catch (IOException ex) {
                Logger.getLogger(Astro_Dots_Tools.class.getName()).log(Level.SEVERE, null, ex);
            }
            IJ.showProgress(n, r.getImageCount());
            stack.addSlice("" + (n + 1), ip);
        }
        ImagePlus imgStack = new ImagePlus(name+"_C "+channel, stack);
        imgStack.setCalibration(cal);
        return imgStack;
    }
    
    /**
     * Dialog ask for channels order and if needed spatial calibration
     * @param channels
     * @param showCal
     * @return ch
     */
    public static ArrayList dialog(String[] channels, boolean showCal, Calibration cal, boolean askParams) {
        ArrayList ch = new ArrayList();
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.addChoice("DAPI : ", channels, channels[0]);
        gd.addChoice("Astrocyte : ", channels, channels[1]);
        gd.addChoice("Dots : ", channels, channels[2]);
        gd.addCheckbox("Do bleaching correction", true);
        if (askParams) {
            gd.addNumericField("Min ratio Intensity", 1, 3);
            gd.addNumericField("Max ratio Intensity", 1, 3);
        }
        if (showCal) {
            gd.addNumericField("XY pixel size : ", cal.pixelWidth, 3);
            gd.addNumericField("Z pixel size : ", 1, 3);
        }
        gd.showDialog();
        ch.add(0, gd.getNextChoiceIndex());
        ch.add(1, gd.getNextChoiceIndex());
        ch.add(2, gd.getNextChoiceIndex());
        doBleachCorr = gd.getNextBoolean();
        if (askParams) {
            minThreshold = gd.getNextNumber();
            maxThreshold = gd.getNextNumber();
        }
        if (showCal) {
            cal.pixelWidth = gd.getNextNumber();
            cal.pixelDepth = gd.getNextNumber();
        }
        return(ch);
    }
    
    /**Gaussian filter 
     * 
     * @param img
     * @param size
     */ 
    public static void gs_filter(ImagePlus img, double size) {
        GaussianBlur gaussian = new GaussianBlur();
        for (int s = 1; s <= img.getNSlices(); s++) {
            img.setZ(s);
            gaussian.blurGaussian(img.getProcessor(), size, size, 0.02);
            img.updateAndDraw();
        }
    }
    
    /*Median filter 
     * 
     * @param img
     * @param size
     */ 
    public static void median_filter(ImagePlus img, double size) {
        RankFilters median = new RankFilters();
        for (int s = 1; s <= img.getNSlices(); s++) {
            img.setZ(s);
            median.rank(img.getProcessor(), size, RankFilters.MEDIAN);
            img.updateAndDraw();
        }
    }
    
    
    /**
     * Find Z with max intensity in stack
     * @param img
     * @return z
     */
    
    private static int find_max(ImagePlus img) {
        double max = 0;
        int zmax = 0;
        for (int z = 1; z <= img.getNSlices(); z++) {
            ImageProcessor ip = img.getStack().getProcessor(z);
            ImageStatistics statistics = new ImageStatistics().getStatistics(ip, ImageStatistics.MEAN, img.getCalibration());
            double meanInt = statistics.mean;
            if (meanInt > max) {
                max = meanInt;
                zmax = z;
            }
        }
        return(zmax);
    }
    
    public static Object3D nucleusSelect(ImagePlus img, Objects3DPopulation nucPop) {
        int nucSelected = find_astroNuc(nucPop, img);
        ImageHandler imgObjNuc = ImageHandler.wrap(img).createSameDimensions();
        ImageHandler imgObjNucSelected = ImageHandler.wrap(img).createSameDimensions();
        ImagePlus imgAstro = doZProjection(img.duplicate(), ZProjector.MAX_METHOD);
        IJ.run(imgAstro, "Enhance Contrast", "saturated=0.35");
        String[] nucIndex = new String[nucPop.getNbObjects()];
        boolean[] values = new boolean[nucPop.getNbObjects()];
        for (int i = 0; i < nucPop.getNbObjects(); i++) {
            Object3D obj = nucPop.getObject(i);
            nucIndex[i] = Integer.toString(i);
            values[i] = false;
            if (i == nucSelected) {
                obj.draw(imgObjNucSelected);
                labelsObject(obj, imgObjNucSelected.getImagePlus(), i-1, 255);
                values[i] = true;
            }
            else {
                obj.draw(imgObjNuc);
                labelsObject(obj, imgObjNuc.getImagePlus(), i-1, 255);
            }
        }
        ImagePlus imhProj = doZProjection(imgObjNuc.getImagePlus(), ZProjector.MAX_METHOD);
        IJ.run(imhProj, "Enhance Contrast", "saturated=0.35");
        imgObjNuc.closeImagePlus();
        ImagePlus imhProjSelected = doZProjection(imgObjNucSelected.getImagePlus(), ZProjector.MAX_METHOD);
        IJ.run(imhProjSelected, "Enhance Contrast", "saturated=0.35");
        
        ImagePlus[] imgColors = {imhProj, imhProjSelected, null, imgAstro, null};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        flush_close(imhProj);
        flush_close(imhProjSelected);
        
        imgObjects.show("Nucleus");
        GenericDialog gd = new GenericDialog("Choose astrocyte nucleus");
        gd.addCheckboxGroup(nucPop.getNbObjects(), 1, nucIndex, values);
        //gd.addRadioButtonGroup("Nucleus", nucIndex, nucPop.getNbObjects(), 0, Integer.toString(nucSelected));
        gd.showDialog();
        Vector<?> checkboxes = gd.getCheckboxes();
        //int nucChoice = Integer.parseInt(gd.getNextRadioButton());
        imgObjects.hide();
        flush_close(imgObjects);
        
        // if more than one object selected merge objects.
        ArrayList<Object3D> objSelected = new ArrayList<>();

        for (int n = 0; n < checkboxes.size(); n++) {
            Checkbox nucCheck = (Checkbox) checkboxes.get(n);
            if (nucCheck.getState()) {
                objSelected.add(nucPop.getObject(n));
                //System.out.println(n);
            }
        }
        Object3DVoxels objVoxel = new Object3DVoxels(imgObjNuc);
        objVoxel.addVoxelsUnion(objSelected);
        return(objVoxel);
    }
    
    
    
    
    /**
     * Threshold images and fill holes
     * @param img
     * @param thMed
     * @param fill 
     */
    public static void threshold(ImagePlus img, AutoThresholder.Method thMed, boolean fill) {
        //  Threshold and binarize
        img.setZ(find_max(img));
        img.updateAndDraw();
        IJ.setAutoThreshold(img, thMed.toString()+" dark");
        Prefs.blackBackground = false;
        IJ.run(img, "Convert to Mask", "method="+thMed.toString()+" background=Dark");
        if (fill) {
            IJ.run(img,"Fill Holes", "stack");
        }
    }
    
    
    public static Objects3DPopulation getPopFromImage(ImagePlus img, Calibration cal) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        labels.setCalibration(cal);
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    }
    
    
    /**
     * Clear out side roi
     * @param img
     * @param roi
     */
    public static void clearOutSide(ImagePlus img, Roi roi) {
        for (int n = 1; n <= img.getNSlices(); n++) {
            ImageProcessor ip = img.getImageStack().getProcessor(n);
            ip.setRoi(roi);
            ip.setBackgroundValue(0);
            ip.fillOutside(roi);
        }
        img.updateAndDraw();
    }
                            
    /**
     * Size filter objects
     * remove touching border
     * 
     * @param min
     * @param max
     * @param objects
     * @param img
     * @param border
    */
    public static void objectsSizeFilter(double min, double max, Objects3DPopulation objects, ImagePlus img, boolean border) {
        ImageHandler imh = ImageInt.wrap(img).createSameDimensions();
        Object3D obj = null;
        boolean remove = false;
        if (objects.getNbObjects() > 0) {
            for (int n = 0; n < objects.getNbObjects(); n++) {
                remove = false;
                obj = objects.getObject(n);
                double vol = obj.getVolumeUnit();
                //System.out.println(vol);
                // remove if touching border
                if (border) {
                    if (obj.touchBorders(imh, false)) {
                        remove = true;
                    }
                }
                // Size filter remove small
                if ((vol < min) || (vol > max)) {
                    remove = true;
                }
                if (remove) {
                    objects.removeObject(n);
                    n--;
                }
            }
        }
    }
    
    
    private static ImagePlus doZProjection(ImagePlus img, int param) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(param);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
    /**
    * Find background image intensity
    * Z project min intensity
    * read mean intensity
    * substract to image
    * @param img
    */
    public static double find_background(ImagePlus img) {
      img.deleteRoi();
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      double baseLineInt = imp.getStatistics().mean;
      flush_close(imgProj);
      System.out.println("Background = "+baseLineInt);
      //IJ.run(img, "Subtract...", "value="+baseLineInt+" stack");
      return(baseLineInt);
    }
    
    /**
     * Bleach correction
     * @param img 
     */
    public static void doBleachCorrection (ImagePlus img) {
        double bg = find_background(img);
        BleachCorrection_SimpleRatio bleachCorr = new BleachCorrection_SimpleRatio(img, bg);
        ImagePlus imgcorrectBleach = bleachCorr.correctBleach();
        img = imgcorrectBleach.duplicate();
        imgcorrectBleach.close();
    }
    
    /**
     * Find roi area and perimeter
     * @param roi
     * @param img
     * @return perimeter, area
     */
    private static double[] roi_parameter(Roi roi, ImagePlus img) {
        double[] param = new double[2];
        ImageProcessor ip = img.getProcessor();
        ip.setRoi(roi);
        ImageStatistics stats = ImageStatistics.getStatistics(ip, 1, img.getCalibration());
        param[0] = roi.getLength();
        param[1] = stats.area;
        return(param);
    }
    
    /**
     * Find astrocyte nucleus keep nucleus that have the higher intensity in astrocyte channel
     * @param nucPop
     * @param imgAstro
     * @return nucAstro
     */
    
    public static int find_astroNuc( Objects3DPopulation nucPop, ImagePlus imgAstro) {
        double maxInt = 0;
        int index = 0;      
        for (int o = 0; o < nucPop.getNbObjects(); o++) {
            Object3D obj = nucPop.getObject(o);
            double astroInt = obj.getPixMeanValue(ImageHandler.wrap(imgAstro));
            //System.out.println(astroInt);
            if (astroInt > maxInt) {
                maxInt = astroInt;
                index = o;
            }
        }
        return(index);
    }
    
    public static ImagePlus watershedSplit(ImagePlus imgBinmask, int rad) {
        int nbCpus = ThreadUtil.getNbCpus();
        float resXY = 1;
        float resZ = 1;
        float radXY = rad;
        float radZ = rad;
        Calibration cal = imgBinmask.getCalibration();
        if (cal != null) {
            resXY = (float) cal.pixelWidth;
            resZ = (float) cal.pixelDepth;
            radZ = radXY * (resXY / resZ);
        }
        IJ.showStatus("Computing EDT");
        ImageInt imgMask = ImageInt.wrap(imgBinmask);
        ImageFloat edt = EDT.run(imgMask, 0, resXY, resZ, false, nbCpus);
        ImageHandler edt16 = edt.convertToShort(true);
        ImagePlus edt16Plus = edt16.getImagePlus();
        GaussianBlur3D.blur(edt16Plus, 2.0, 2.0, 2.0);
        edt16 = ImageInt.wrap(edt16Plus);
        edt16.intersectMask(imgMask);
        // seeds
        ImageHandler seedsImg;
        seedsImg = FastFilters3D.filterImage(edt16, FastFilters3D.MAXLOCAL, radXY, radXY, radZ, 0, false);
        
        IJ.showStatus("Computing watershed");
        Watershed3D water = new Watershed3D(edt16, seedsImg, 0, 0);
        ImagePlus imp = water.getWatershedImage3D().getImagePlus();
        IJ.setThreshold(imp, 1, 65535);
        Prefs.blackBackground = false;
        IJ.run(imp, "Convert to Mask", "method=Default background=Dark");
        return(imp);
    }
    
    /**
     * Find dots population
     * @param imgDots
     * @return dotsPop
     */
    
    public static Objects3DPopulation find_dots(ImagePlus imgDots, Roi roi) {
        ImagePlus imp = imgDots.duplicate();
        median_filter(imgDots, 1.5);
        IJ.run(imp, "Difference of Gaussians", "  sigma1=3 sigma2=1 enhance stack");
        threshold(imp, AutoThresholder.Method.Li, false);
        if (roi != null) {
            imp.setRoi(roi);
            IJ.run(imp, "Clear Outside", "stack");
        }
        Objects3DPopulation dotsPop = getPopFromImage(imp, imgDots.getCalibration());
        return(dotsPop);
    }
    
    
    // write object labels
    public static void labelsObject (Object3D obj, ImagePlus img, int number, int color) {
        Font tagFont = new Font("SansSerif", Font.PLAIN, 24);
        int[] box = obj.getBoundingBox();
        int z = (int)obj.getCenterZ();
        int x = box[0] - 2;
        int y = box[2] - 2;
        img.setSlice(z+1);
        ImageProcessor ip = img.getProcessor();
        ip.setFont(tagFont);
        ip.setColor(color);
        ip.drawString(Integer.toString(number+1), x, y);
        img.updateAndDraw();    
    }
    
   /**
    * draw population
    * nucleus in blue
    * over baseline value = 1 in green
    * under baseline value = 0 in red
    * @param nucObj
    * @param dotsPop
    * @param img
     * @param imgDir
     * @param imgName
     * @param roi
    */
    public static void tagsObjects(Object3D nucObj, Objects3DPopulation dotsPop, ImagePlus img, String imgDir, String imgName, int roi) {  
        ImageHandler imgObjNuc = ImageHandler.wrap(img).createSameDimensions();
        nucObj.draw(imgObjNuc, 255);                        
        ImageHandler imgObjDotsNotAstro = ImageHandler.wrap(img).createSameDimensions();
        ImageHandler imgObjDotsFineProcess = ImageHandler.wrap(img).createSameDimensions();
        ImageHandler imgObjDotsLargeProcess = ImageHandler.wrap(img).createSameDimensions();
        for (int n = 0; n < dotsPop.getNbObjects(); n++) {
            Object3D obj = dotsPop.getObject(n);
            if (obj.getValue() == 0)
                obj.draw(imgObjDotsNotAstro, 255);
            else if (obj.getValue() == 1)
                obj.draw(imgObjDotsFineProcess, 255);
            else
                obj.draw(imgObjDotsLargeProcess, 255);
        }
        // save image for objects population
        ImagePlus[] imgColors = {imgObjDotsNotAstro.getImagePlus(), imgObjDotsLargeProcess.getImagePlus(), imgObjNuc.getImagePlus(),
            img,null,imgObjDotsFineProcess.getImagePlus()};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(img.getCalibration());
        IJ.run(imgObjects, "Enhance Contrast", "saturated=0.35");
//        imgObjects.updateAndDraw();
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(imgDir + imgName + "_Astro-"+(roi+1)+"-Objects.tif"); 
        flush_close(imgObjects);
    }
    
    
   
// Flush and close images
    public static void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    } 
         
    /**
    * Compute dots parameters
    * @param roi roi of astrocyte
    * @param nucObj nucleus
    * @param dotsPop dots population
     * @param imgDots
    * @param imgAstro read dots intensity
     * @param results buffer
     * @throws java.io.IOException 
    **/
    public static void compute_parameters(Roi roi, Object3D nucObj, Objects3DPopulation dotsPop, ImagePlus imgDots, ImagePlus imgAstro, 
            BufferedWriter results) throws IOException {
        IJ.showStatus("Computing nucleus parameters ....");

        // measure roi volume
        double[] params = roi_parameter(roi, imgAstro);
        double perimeter = params[0];
        double area = params[1];
        //System.out.println(area+", "+perimeter);
        imgAstro.deleteRoi();
        double astroVol = area * imgAstro.getNSlices() * imgAstro.getCalibration().pixelDepth;
        double circularity = 4*Math.PI*(area/Math.pow(perimeter,2));
        
        // measure dots integrated intensity and volume
        double bgInt = find_background(imgAstro);
        double minThresholdInt =  bgInt*minThreshold;
        double maxThresholdInt =  bgInt*maxThreshold;
        for (int n = 0; n < dotsPop.getNbObjects(); n++) {
            Object3D dotObj = dotsPop.getObject(n);
            double volDot = dotObj.getVolumeUnit();
            double intDotAstroimg = dotObj.getIntegratedDensity(ImageHandler.wrap(imgAstro));
            /* classify dots
            * dots 0  dots < minThresholdInt and dist > 2 (not in astro (red))
            * dots 1  minThresholdInt =< dots < maxThresholdInt and dist > 2 (fine processes (magenta))
            * dots 2 dots >= maxThresholdInt or dist < 2 (larger processes (green))
            */
            double dist = 0;
            if (nucObj.includes(dotObj))
                dist = 0;
            else
                dist = dotObj.distBorderUnit(nucObj);
            // dots 0
            if ((intDotAstroimg < minThresholdInt) && (dist < 2))
                dotObj.setValue(0);
            // dost 1
            else if ((intDotAstroimg >= minThresholdInt) && (intDotAstroimg < maxThresholdInt) && (dist > 2))
                dotObj.setValue(1);
            // dots 2
            else
                dotObj.setValue(2);
            
            // Compute nucleus and soma (enlarge nucleus 2µm ) volume
            double nucVol = nucObj.getVolumeUnit();
            double somaVol = 4/3*Math.PI*(3*Math.pow(nucObj.getFeret(),2)+3*nucObj.getFeret()+2);
            results.write(roi.getName() + "\t" + astroVol+"\t"+circularity+"\t"+nucVol+"\t"+somaVol+"\t"+n+"\t"+volDot+"\t"+intDotAstroimg+"\t"+dist+"\t"+dotObj.getValue()+"\t"+bgInt+"\n");
            results.flush();
        }
    }
    
    /**
     * Compute and write global image parameters
     */
    public static void compute_Image_parameters(String outDir, String imageName) throws IOException {
        File directory = new File(outDir);
        FileFilter filter = new RegexFileFilter(imageName+"_Astro_.*results.xls");
        File[] files = directory.listFiles(filter);
        FileWriter fileResults = new FileWriter(outDir + imageName + "_results.xls");
        BufferedWriter outPutResults = new BufferedWriter(fileResults);
        /* write headers
        * largeBrDots (dots in large branches) = dots 2
        * fineBrDots (dots in fine branches) = dots 1
        * dotsNoinAstro (dots outiside astrocyte) = dots 0
        * totalDots (all dots) = dots0 + dots1 + dots2
        * dotsinAstro (dots in astrocyte) = largeBrDots + fineBrDots
        * dotsinSoma (dots inside soma) = dots with dist to nucleus < 2
        * dotsinNuc (dots in nucleus) = dost with dist = 0
        * dotsdensinAstro (dots density in astrocyte) = dotsinAstro / astroVolume
        * percDotsNotinAstro (% dots not in astrocyte = dotsNotinAstro / totalDots
        * percDotsinSoma (% dots in soma) = dotsinSoma/totalDotsAstro
        * perDostFineBr (% dost in fine branches) = fineBrDots / totalDotsAstro
        */
        outPutResults.write("Roi Name\tAstrocyte Volume\tAstrocyte circularity\tDensity dots in Astro"
                + "\t% of dots not in astro\t% of dots in soma\t% of dots in fine processes\tDots mean distance\tSD distance"
                + "\tDots mean intensity in Astro\tSD intensity in astro\n");
        outPutResults.flush();
        Arrays.sort(files);
    for (File file : files) {
        double largeBrDots = 0, fineBrDots = 0, dotsNoinAstro = 0, totalDots, dotsinAstro, dotsinSoma = 0, dotsinNuc = 0, meanDotsinAstroInt,  sdDotsinAstroInt;
        double astroVol, nucVol, somaVol, astroCirc, dotsdensinAstro, percDotsNotinAstro, percDotsinSoma, perDostFineBr;

        ResultsTable rt = new ResultsTable().open(file.toString());
        String roiName =  rt.getStringValue("Roi Name", 0);
        astroVol = rt.getValue("Astrocyte Volume", 0);
        astroCirc = rt.getValue("Astrocyte circularity", 0);
        nucVol = rt.getValue("Nucleus Vol",0);
        somaVol = rt.getValue("Soma Vol", 0);
        DescriptiveStatistics statsDist = new DescriptiveStatistics();
        DescriptiveStatistics statsDotsinAstroInt = new DescriptiveStatistics();
        totalDots = rt.getCounter();
        
        for (int r = 0; r < totalDots; r++) {
            double nucDist = rt.getValue("Distance to nucleus border", r);
            int Dottype = (int)rt.getValue("Dots type", r);
            switch (Dottype) {
                case 0 :
                    if (nucDist > 2)
                        dotsNoinAstro++;
                    else
                        largeBrDots++;
                    break;
                case 1 :
                    if (nucDist > 2)
                        fineBrDots++;
                    break;
                case 2 :    
                    largeBrDots++;
                    break;
            }        
            if (nucDist < 2)
                dotsinSoma++;
            if (nucDist == 0)
                dotsinNuc++;
            statsDist.addValue(nucDist);
            double dotsinAstroInt = rt.getValue("Dot Norm intDensity in Astro", r);
            statsDotsinAstroInt.addValue(dotsinAstroInt);
        }
        dotsinAstro = largeBrDots + fineBrDots;
        dotsdensinAstro = dotsinAstro / astroVol;
        percDotsNotinAstro = 100*(dotsNoinAstro / totalDots);
        percDotsinSoma = 100*(dotsinSoma / dotsinAstro);
        perDostFineBr = 100*(fineBrDots / dotsinAstro);
        meanDotsinAstroInt = statsDotsinAstroInt.getMean();
        sdDotsinAstroInt = statsDotsinAstroInt.getStandardDeviation();
        
        outPutResults.write(roiName+ "\t" + astroVol + "\t" + astroCirc + "\t" + dotsdensinAstro + "\t" + percDotsNotinAstro + "\t" + percDotsinSoma
                + "\t" + perDostFineBr + "\t" + statsDist.getMean() + "\t" + statsDist.getStandardDeviation() + "\t" + meanDotsinAstroInt + "\t" + sdDotsinAstroInt +"\n");
        outPutResults.flush();
        rt.reset();
    }
        fileResults.close();
    }
}
