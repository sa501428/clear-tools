package cli.utils.sift;

import javastraw.reader.block.ContactRecord;

import java.util.List;

public class BoundingBoxWithContacts {

    private final static int buffer = 10;
    private final List<ContactRecord> contacts;
    private int minR, maxR, minC, maxC;

    public BoundingBoxWithContacts(List<ContactRecord> contacts) {
        this.contacts = contacts;
        minR = contacts.get(0).getBinX();
        minC = contacts.get(0).getBinY();
        maxR = minR;
        maxC = minC;
        setBounds();
    }

    private void setBounds() {
        for (ContactRecord contact : contacts) {
            minR = Math.min(minR, contact.getBinX() - buffer);
            minC = Math.min(minC, contact.getBinY() - buffer);

            maxR = Math.max(maxR, contact.getBinX() + buffer);
            maxC = Math.max(maxC, contact.getBinY() + buffer);
        }
        if (minR < 0) minR = 0;
        if (minC < 0) minC = 0;
    }
}
