package cli.utils.data;

import javastraw.expected.Welford;
import javastraw.expected.Zscore;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

public class InterLogZscoreExpectedModel {

    private final Welford welford = new Welford();

    public InterLogZscoreExpectedModel(MatrixZoomData zd, NormalizationType norm) {
        Iterator<ContactRecord> records = zd.getNormalizedIterator(norm);
        while (records.hasNext()) {
            ContactRecord record = records.next();
            if (record.getCounts() > 0) {
                welford.addValue(Math.log(record.getCounts()));
            }
        }
    }


    public List<ContactRecord> getAllContactsZscore(MatrixZoomData zd, NormalizationType norm, float threshold) {
        Iterator<ContactRecord> records = zd.getNormalizedIterator(norm);
        Zscore zscore = welford.getZscore();
        List<ContactRecord> topRecords = new LinkedList<>();
        while (records.hasNext()) {
            ContactRecord record = records.next();
            if (record.getCounts() > 0) {
                double z = zscore.getZscore(Math.log(record.getCounts()));
                if (z > threshold) {
                    topRecords.add(record);
                }
            }
        }
        return topRecords;
    }
}
