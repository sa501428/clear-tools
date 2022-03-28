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

import cli.Main;
import jargs.gnu.CmdLineParser;

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
    private final Option normalizationTypeOption = addStringOption('k', "norm");
    private final Option cutoffOption = addIntegerOption("cutoff");

    public CommandLineParser() {
    }

    /*
     * convert Options to Objects or Primitives
     */

    private boolean optionToBoolean(Option option) {
        Object opt = getOptionValue(option);
        return opt != null && (Boolean) opt;
    }

    private int optionToInteger(Option option, int defaultValue) {
        Object opt = getOptionValue(option);
        return opt == null ?  defaultValue : ((Number) opt).intValue();
    }

    private String optionToString(Option option) {
        Object opt = getOptionValue(option);
        return opt == null ? Main.DEFAULT_NORMALIZATION : opt.toString();
    }

    /*
     * Actual parameters
     */

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

    public int getCutoffOption(){
        return optionToInteger(cutoffOption, Main.DEFAULT_CUTOFF);
    }

    public int getResolutionOption(){
        return optionToInteger(resolutionsOption, Main.DEFAULT_RESOLUTION);
    }

    public String getNormalizationStringOption() {
        return optionToString(normalizationTypeOption);
    }


}