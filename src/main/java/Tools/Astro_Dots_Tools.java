package Tools;

import static Astrocytes.Astro_Dots.meanDotsSize;
import static Astrocytes.Astro_Dots.stdMeanDotsSize;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.*;
import ij.plugin.Duplicator;
import ij.plugin.GaussianBlur3D;
import ij.plugin.RGBStackMerge;
import ij.plugin.Thresholder;
import ij.plugin.ZProjector;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.RankFilters;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.awt.Checkbox;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
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
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import sc.fiji.localThickness.Clean_Up_Local_Thickness;
import sc.fiji.localThickness.EDT_S1D;
import sc.fiji.localThickness.Local_Thickness_Parallel;


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
public static boolean useDeconvImg = false;
public static String thMethod = "";
public static double bg = 0;
public static double stdBg = 0;
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
        String[] thMethods = new Thresholder().methods;
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.addChoice("DAPI : ", channels, channels[0]);
        gd.addChoice("Astrocyte : ", channels, channels[1]);
        gd.addChoice("Dots : ", channels, channels[2]);
        if (askParams) {
            gd.addNumericField("Mean dot Volume", meanDotsSize, 3);
            gd.addNumericField("Std dot Volume", stdMeanDotsSize, 3);
        }
        if (showCal) {
            gd.addNumericField("XY pixel size : ", cal.pixelWidth, 3);
            gd.addNumericField("Z pixel size : ", 1, 3);
        }
        gd.addChoice("Threshold Methods", thMethods, thMethods[5]);
        gd.addCheckbox("Use deconvolved image ", false);
        gd.showDialog();
        ch.add(0, gd.getNextChoiceIndex());
        ch.add(1, gd.getNextChoiceIndex());
        ch.add(2, gd.getNextChoiceIndex());
        thMethod = thMethods[gd.getNextChoiceIndex()];
        if (askParams) {
            meanDotsSize = gd.getNextNumber();
            stdMeanDotsSize = gd.getNextNumber();
        }
        if (showCal) {
            cal.pixelWidth = gd.getNextNumber();
            cal.pixelDepth = gd.getNextNumber();
        }
        useDeconvImg = gd.getNextBoolean();
        if(gd.wasCanceled())
            ch = null;
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
                labelsObject(obj, imgObjNucSelected.getImagePlus(), i, 255, 24);
                values[i] = true;
            }
            else {
                obj.draw(imgObjNuc);
                labelsObject(obj, imgObjNuc.getImagePlus(), i, 255, 24);
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
    public static void threshold(ImagePlus img, String thMed, boolean fill) {
        //  Threshold and binarize
        img.setZ(find_max(img));
        img.updateAndDraw();
        IJ.setAutoThreshold(img, thMed + " dark");
        Prefs.blackBackground = false;
        IJ.run(img, "Convert to Mask", "method="+thMed+" background=Dark");
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
    public static void find_background(ImagePlus img) {
      img.deleteRoi();
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      bg = imp.getStatistics().mean;
      stdBg = imp.getStatistics().stdDev;
      flush_close(imgProj);
    }
    
    
    /**
     * Find roi volume
     * @param roi
     * @param img
     * @return volume
     */
    private static double roi_volume(Roi roi, ImagePlus imgAstro) {
        ImagePlus img = new Duplicator().run(imgAstro, 1, 1);
        ImageProcessor ip = img.getProcessor();
        ip.setRoi(roi);
        ImageStatistics stats = ImageStatistics.getStatistics(ip, Measurements.MEAN, img.getCalibration());
        double volume = stats.area * img.getCalibration().pixelDepth * img.getNSlices();
        flush_close(img);
        return(volume);
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
    
    /**
     * compute local thickness
     * @param imgAstro
     * @return astroMap
    **/
    public static ImagePlus localThickness3D (ImagePlus imgAstro) {
        ImagePlus img = imgAstro.duplicate();
        median_filter(img, 1);
        threshold(img, "Li", false);
        EDT_S1D edt = new EDT_S1D();
        edt.runSilent = true;
        edt.thresh = 1;
        edt.inverse = false;
        edt.showOptions = false;
        edt.setup("", img);
        edt.run(img.getProcessor());
        Local_Thickness_Parallel locThk = new Local_Thickness_Parallel();
        locThk.runSilent = true;
        locThk.setup("", edt.getResultImage());
        locThk.run(edt.getResultImage().getProcessor());
        Clean_Up_Local_Thickness cleanUp = new Clean_Up_Local_Thickness();
        cleanUp.runSilent = true;
        cleanUp.setup("", locThk.getResultImage());
        cleanUp.run(locThk.getResultImage().getProcessor());
        
        ImagePlus astroMap = cleanUp.getResultImage();
//        IJ.run(img, "Options...", "iterations=1 count=1 do=Erode stack");
//        String imageName = img.getTitle()+"_LocThk";
//        IJ.showStatus("Computing local thickness ...");
//        //IJ.run(img, "Local Thickness (masked, calibrated, silent)", "");
//        IJ.run(img, "Local Thickness (complete process)", "threshold=1");
//        while (WindowManager.getImage(imageName) == null)
//            IJ.wait(10);
//        ImagePlus astroMap = WindowManager.getImage(imageName);
//        astroMap.hide();
        astroMap.setCalibration(img.getCalibration());
        flush_close(img);
        return(astroMap);
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
        ImageFloat edt = EDT.run(imgMask, 1, resXY, resZ, false, nbCpus);
        ImageHandler edt16 = edt.convertToShort(true);
        ImagePlus edt16Plus = edt16.getImagePlus();
        GaussianBlur3D.blur(edt16Plus, 6.0, 6.0, 6.0);
        edt16 = ImageInt.wrap(edt16Plus);
        edt16.intersectMask(imgMask);
        // seeds
        ImageHandler seedsImg;
        seedsImg = FastFilters3D.filterImage(edt16, FastFilters3D.MAXLOCAL, radXY, radXY, radZ, 0, false);
        IJ.showStatus("Computing watershed");
        Watershed3D water = new Watershed3D(edt16, seedsImg, 0, 1);
        ImagePlus imp = water.getWatershedImage3D().getImagePlus();
        IJ.setThreshold(imp, 1, 65535);
        Prefs.blackBackground = false;
        IJ.run(imp, "Convert to Mask", "method=Default background=Dark");
        flush_close(imgMask.getImagePlus());
        flush_close(edt.getImagePlus());
        flush_close(edt16.getImagePlus());
        flush_close(edt16Plus);
        flush_close(seedsImg.getImagePlus());
        return(imp);
    }
    
    /**
     * Find dots population
     * @param imgDots
     * @param roi
     * @param deconv
     * @return dotsPop
     */
    
    public static Objects3DPopulation find_dots(ImagePlus imgDots, Roi roi, boolean deconv) {
        ImagePlus imp = imgDots.duplicate();
        median_filter(imgDots, 1);
        if (!deconv)
            IJ.run(imp, "Difference of Gaussians", "  sigma1=3 sigma2=1 enhance stack");
        threshold(imp, thMethod, false);   
        if (roi != null) {
            imp.setRoi(roi);
            IJ.run(imp, "Clear Outside", "stack");
        }
        Objects3DPopulation dotsPop = getPopFromImage(imp, imgDots.getCalibration());
        return(dotsPop);
    }
    
    
    // write object labels
    public static void labelsObject (Object3D obj, ImagePlus img, int number, int color, int size) {
        Font tagFont = new Font("SansSerif", Font.PLAIN, size);
        int[] box = obj.getBoundingBox();
        int z = (int)obj.getCenterZ();
        int x = box[0] - 2;
        int y = box[2] - 2;
        img.setSlice(z+1);
        ImageProcessor ip = img.getProcessor();
        ip.setFont(tagFont);
        ip.setColor(color);
        ip.drawString(Integer.toString(number), x, y);
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
        ImagePlus imAstro = img.duplicate();
        ImageHandler imgObjNuc = ImageHandler.wrap(imAstro).createSameDimensions();
        nucObj.draw(imgObjNuc, 255);                        
        ImageHandler imgObjDotsNotAstro = ImageHandler.wrap(imAstro).createSameDimensions();
        ImageHandler imgObjDotsFineProcess = ImageHandler.wrap(imAstro).createSameDimensions();
        ImageHandler imgObjDotsLargeProcess = ImageHandler.wrap(imAstro).createSameDimensions();
        for (int n = 0; n < dotsPop.getNbObjects(); n++) {
            Object3D obj = dotsPop.getObject(n);
            switch (obj.getValue()) {
                case 0:
                    obj.draw(imgObjDotsNotAstro, 255);
                    break;
                case 1:
                    obj.draw(imgObjDotsFineProcess, 255);
                    break;
                default:
                    obj.draw(imgObjDotsLargeProcess, 255);
                    break;
            }
        }
        // save image for objects population
        ImagePlus[] imgColors = {imgObjDotsNotAstro.getImagePlus(), imgObjDotsLargeProcess.getImagePlus(), imgObjNuc.getImagePlus(),
            imAstro, imgObjDotsFineProcess.getImagePlus()};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(img.getCalibration());
        IJ.run(imgObjects, "Enhance Contrast", "saturated=0.35");
        imgObjects.updateAndDraw();
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
    * classify dots
    * @param nucObj nucleus
    * @param dotsPop dots population
    * @param imgAstro read dots intensity
     * @param astroMap
    **/
    public static void classify_dots(Object3D nucObj, Objects3DPopulation dotsPop, ImagePlus imgAstro, ImagePlus astroMap) {
        IJ.showStatus("Classify dots ....");
        find_background(imgAstro);
        double BgThresholdInt = bg + stdBg;
        // measure dots mean intensity and volume
        for (int n = 0; n < dotsPop.getNbObjects(); n++) {
            Object3D dotObj = dotsPop.getObject(n);
            double astroDiameter = dotObj.getPixMaxValue(ImageHandler.wrap(astroMap)) * imgAstro.getCalibration().pixelWidth;
            double meanIntDotAstroimg = dotObj.getPixMeanValue(ImageHandler.wrap(imgAstro));
            /* classify dots
            * dots 0  dots <= minThresholdInt and dist > 2 (not in astro (red))
            * dots 1  dots > minThresholdInt  and astroDiameter <= spatial res in Z and dist > 2 (fine processes (magenta))
            * dots 2 dots (larger processes (green))
            */
            double distNuc = 0;
            if (nucObj.includes(dotObj))
                distNuc = 0;
            else
                distNuc = dotObj.distBorderUnit(nucObj);
            // dots 0
            if ((meanIntDotAstroimg <= BgThresholdInt) && (distNuc > 2))
                dotObj.setValue(0);
            // dost 1
            
            else if ((meanIntDotAstroimg > BgThresholdInt) && (astroDiameter <= imgAstro.getCalibration().pixelDepth * 2) 
                    && (distNuc > 2))
                dotObj.setValue(1);
            // dots 2
            else
                dotObj.setValue(2);
        }
    }
    
    /**
     * Compute and write global image parameters
     * @param roi
     * @param imgAstro
     * @param imgAstroMap
     * @param nucObj
     * @param dotsPop
     * @param outDir
     * @param imageName
     * @throws java.io.IOException
     */
    public static void compute_Image_parameters(Roi roi, ImagePlus imgAstro, ImagePlus imgAstroMap, Object3D nucObj, Objects3DPopulation dotsPop,
            BufferedWriter results, String imageName) throws IOException {
        /* write headers
        * largeBrDots (dots in large branches) = dots 2
        * fineBrDots (dots in fine branches) = dots 1
        * dotsNoinAstro (dots outiside astrocyte) = dots 0
        * totalDots (all dots) = dots0 + dots1 + dots2
        * dotsinAstro (dots in astrocyte) = largeBrDots + fineBrDots
        * dotsinSoma (dots inside soma) = dots with dist to nucleus < 2
        * dotsinNuc (dots in nucleus) = dots with dist = 0
        * dotsdensinAstro (dots density in astrocyte) = dotsinAstro / astroVolume
        * percDotsNotinAstro (% dots not in astrocyte = dotsNotinAstro / totalDots
        * percDotsinSoma (% dots in soma) = dotsinSoma/totalDotsAstro
        * perDostFineBr (% dost in fine branches) = fineBrDots / totalDotsAstro
        */
        
        // measure roi volume
        double astroVol = roi_volume(roi, imgAstro);
        // Compute nucleus volume
        // add to astro parameters
        double nucVol = nucObj.getVolumeUnit(); 
        float dotsNoinAstro = 0;
        float dotsinFineBr = 0;
        float dotsinLargeBr = 0;
        float dotsinSoma= 0;
        DescriptiveStatistics dotVolStats = new DescriptiveStatistics();
        DescriptiveStatistics dotDiameterStats = new DescriptiveStatistics();
        DescriptiveStatistics dotMeanIntStats = new DescriptiveStatistics();
        DescriptiveStatistics dotNucDistanceStats = new DescriptiveStatistics();
        String roiName = roi.getName();
        for (int n = 0; n < dotsPop.getNbObjects(); n++) {
            Object3D dotObj = dotsPop.getObject(n);
            double dotVol = dotObj.getVolumeUnit();
            dotVolStats.addValue(dotVol);
            double distNuc = 0;
            if (nucObj.includes(dotObj))
                distNuc = 0;
            else {
                distNuc = dotObj.distCenterBorderUnit(nucObj);
                dotNucDistanceStats.addValue(distNuc);
            }
            int Dottype = dotObj.getValue();
            double astroDiameter = dotObj.getPixMaxValue(ImageHandler.wrap(imgAstroMap));
            if (astroDiameter != 0)
                dotDiameterStats.addValue(astroDiameter);
            double meanDotSizeStd = meanDotsSize+stdMeanDotsSize;
            switch (Dottype) {
                case 0 :
                    if (dotVol >= meanDotSizeStd) {
                        double ratioVol = Math.round(dotVol/meanDotSizeStd);
                        dotsNoinAstro+=ratioVol;
                    }
                    else
                        dotsNoinAstro++;
                    break;
                case 1 :
                    if (dotVol >= meanDotSizeStd) {
                        double ratioVol = Math.round(dotVol/meanDotSizeStd);
                        dotsinFineBr+=ratioVol;
                    }
                    else
                        dotsinFineBr++;
                    break;
                case 2 : 
                    if (dotVol >= meanDotSizeStd) {
                        double ratioVol = Math.round(dotVol/meanDotSizeStd);
                        dotsinLargeBr+= ratioVol;
                    }
                    else {
                        dotsinLargeBr++;
                        break;
                    }
            }
            if (distNuc < 2) 
                if (dotVol >= meanDotSizeStd) {
                    double ratioVol = Math.round(dotVol/meanDotSizeStd);
                    dotsinSoma+=ratioVol;  
                }   
                else 
                    dotsinSoma++;
            
            double IntdotsinAstro = dotObj.getPixMeanValue(ImageHandler.wrap(imgAstro));
            dotMeanIntStats.addValue(IntdotsinAstro);
        }
        float totalDots = dotsinLargeBr + dotsinFineBr + dotsNoinAstro;
        float dotsinAstro = dotsinLargeBr + dotsinFineBr;
        double dotsdensinAstro = dotsinAstro / astroVol;
        double percDotsNotinAstro = 100*(dotsNoinAstro / totalDots);
        double percDotsinSoma = 100*(dotsinSoma / dotsinAstro);
        double perDotsFineBr = 100*(dotsinFineBr / dotsinAstro);
        double perDotsLargeBr = 100*(dotsinLargeBr / dotsinAstro);
        double meanIntDotsinAstro = dotMeanIntStats.getMean();
        double sdIntDotsinAstro = dotMeanIntStats.getStandardDeviation();
        double meanAstroDiameter = dotDiameterStats.getMean();
        double stdAstroDiameter = dotDiameterStats.getStandardDeviation();
        double medAstroDiameter = dotDiameterStats.getPercentile(50); 
        double meanNucDist = dotNucDistanceStats.getMean();
        double stdNucDist = dotNucDistanceStats.getStandardDeviation();
        results.write(imageName + "\t" + roiName + "\t" + bg + "\t" + stdBg + "\t" + astroVol + "\t" + dotsdensinAstro +
                "\t" + percDotsNotinAstro + "\t" + percDotsinSoma + "\t" + perDotsFineBr + "\t" + perDotsLargeBr +
                "\t" + meanIntDotsinAstro + "\t" + sdIntDotsinAstro + "\t"+meanAstroDiameter+"\t"+stdAstroDiameter+
                "\t"+medAstroDiameter+"\t"+meanNucDist+"\t"+stdNucDist+"\n");
        results.flush();
    }
    
}
