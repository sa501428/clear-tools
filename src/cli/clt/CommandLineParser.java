/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2019-2021 Rice University, Baylor College of Medicine, Aiden Lab
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

package cli.clt;

import jargs.gnu.CmdLineParser;
import javastraw.reader.Dataset;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;

/**
 * Command Line Parser for EMT commands
 *
 * @author Muhammad Shamim
 */
public class CommandLineParser extends CmdLineParser {

    private final Option verboseOption = addBooleanOption('v', "verbose");
    private final Option npyOption = addBooleanOption("npy");
    private final Option helpOption = addBooleanOption('h', "help");
    private final Option versionOption = addBooleanOption('V', "version");
    private final Option logOption = addBooleanOption("log");
    private final Option resolutionsOption = addIntegerOption('r', "res");
    private final Option lowResolutionsOption = addIntegerOption("low-res");
    private final Option seedOption = addIntegerOption("seed");
    private final Option normalizationTypeOption = addStringOption('k', "norm");
    private final Option cutoffOption = addIntegerOption("cutoff");
    private final Option percentileOption = addIntegerOption("percentile");
    private final Option windowOption = addIntegerOption('w', "window");
    private final Option minDisValOption = addIntegerOption("min-dist");
    private final Option maxDistValOption = addIntegerOption("max-dist");
    private final Option interChromosomalOption = addBooleanOption("include-inter");
    private final Option onlyOneOption = addBooleanOption("only-one");
    private final Option noOrientationOption = addBooleanOption("no-orientation");
    private final Option aggregateNormalization = addBooleanOption("ag-norm");
    private final Option isLoopAnalysis = addBooleanOption("loop");
    private final Option thresholdOption = addDoubleOption("threshold");
    private final Option attributeOption = addStringOption("attributes");
    private final Option numThreadsOption = addIntegerOption("threads");

    private boolean optionToBoolean(Option option) {
        Object opt = getOptionValue(option);
        return opt != null && (Boolean) opt;
    }

    private int optionToInteger(Option option, int defaultValue) {
        Object opt = getOptionValue(option);
        return opt == null ? defaultValue : ((Number) opt).intValue();
    }

    private double optionToDouble(Option option, double defaultValue) {
        Object opt = getOptionValue(option);
        return opt == null ? defaultValue : ((Number) opt).doubleValue();
    }

    private String optionToString(Option option) {
        Object opt = getOptionValue(option);
        return opt == null ? null : opt.toString();
    }

    private String[] optionToStringArray(Option option) {
        Object opt = getOptionValue(option);
        return opt == null ? null : opt.toString().split(",");
    }

    public boolean getHelpOption() {
        return optionToBoolean(helpOption);
    }

    public boolean getVerboseOption() {
        return optionToBoolean(verboseOption);
    }

    public boolean getNpyOption() {
        return optionToBoolean(npyOption);
    }

    public boolean getVersionOption() {
        return optionToBoolean(versionOption);
    }

    public boolean getLogOption() {
        return optionToBoolean(logOption);
    }

    public int getCutoffOption(int defaultValue) {
        return optionToInteger(cutoffOption, defaultValue);
    }

    public double getThresholdOption(double defaultVal) {
        return optionToDouble(thresholdOption, defaultVal);
    }

    public int getResolutionOption(int defaultVal) {
        return optionToInteger(resolutionsOption, defaultVal);
    }

    public int getLowResolutionOption(int defaultVal) {
        return optionToInteger(lowResolutionsOption, defaultVal);
    }

    public long getSeedOption(int defaultVal) {
        return optionToInteger(seedOption, defaultVal);
    }

    public String getNormalizationStringOption() {
        return optionToString(normalizationTypeOption);
    }

    public boolean getAggregateNormalization() {
        return optionToBoolean(aggregateNormalization);
    }

    public boolean getIsLoopAnalysis() {
        return optionToBoolean(isLoopAnalysis);
    }

    public int getMinDistVal(int val) {
        return optionToInteger(minDisValOption, val);
    }

    public int getMaxDistVal(int val) {
        return optionToInteger(maxDistValOption, val);
    }

    public int getWindowSizeOption(int val) {
        return optionToInteger(windowOption, val);
    }

    public boolean getIncludeInterChromosomal() {
        return optionToBoolean(interChromosomalOption);
    }

    public int getPercentileOption(int val) {
        return optionToInteger(percentileOption, val);
    }

    public NormalizationType getNormOrDefaultScale(Dataset ds) {
        String normType = getNormalizationStringOption();
        if (normType != null && normType.length() > 1) {
            return ds.getNormalizationHandler().getNormTypeFromString(normType);
        } else {
            return NormalizationHandler.SCALE;
        }
    }

    public boolean getOnlyOneOption() {
        return optionToBoolean(onlyOneOption);
    }

    public boolean getNoOrientationFlag() {
        return optionToBoolean(noOrientationOption);
    }

    public String[] getAttributesOption() {
        return optionToStringArray(attributeOption);
    }

    public int getNumThreads(int numThreads) {
        return optionToInteger(numThreadsOption, numThreads);
    }
}