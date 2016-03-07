package Console;

/*Step 2
 *Stat. the distribution of each aa
 */

import java.io.FileInputStream;
import java.io.IOException;

public class Sequence_single_distribution {
    public static void main(String[] args) throws IOException {

        String infilename =
//            "D:\\My Documents\\GThesis\\DatabaseBAK\\std_aa.xml";
        	"D:\\My Documents\\GThesis\\DatabaseBAK\\std_aa_1M_01.xml";

        byte data[] = new byte[1];
        char line[] = new char[120];
        int i = 0, j = 0;
        long stat[] = new long[20];
        float distribution[] = new float[20];
        int k = 0;
        long seq_num = 0;
        long total = 0;

        //Set to Zero
        for (i = 0; i < 20; i++) {
            stat[i] = 0;
        }

        i = 0;

        try {
            FileInputStream infile = new FileInputStream(infilename);

            System.out.print("Reading data...\n");

            while ((infile.read(data, 0, 1)) != -1) {
                if (i > 14 && i < 135 && data[0] > 64 && data[0] < 90) {
                    line[j] = (char) data[0];
                    j++;
                } else if (i == 138 && data[0] == 13 && j == 120) {

                    for (k = 0; k < 120; k++) {
                        switch (line[k]) {
                        case 'G':
                            stat[0]++;
                            break;
                        case 'A':
                            stat[1]++;
                            break;
                        case 'P':
                            stat[2]++;
                            break;
                        case 'V':
                            stat[3]++;
                            break;
                        case 'L':
                            stat[4]++;
                            break;
                        case 'I':
                            stat[5]++;
                            break;
                        case 'M':
                            stat[6]++;
                            break;
                        case 'F':
                            stat[7]++;
                            break;
                        case 'Y':
                            stat[8]++;
                            break;
                        case 'W':
                            stat[9]++;
                            break;
                        case 'S':
                            stat[10]++;
                            break;
                        case 'T':
                            stat[11]++;
                            break;
                        case 'C':
                            stat[12]++;
                            break;
                        case 'N':
                            stat[13]++;
                            break;
                        case 'Q':
                            stat[14]++;
                            break;
                        case 'K':
                            stat[15]++;
                            break;
                        case 'H':
                            stat[16]++;
                            break;
                        case 'R':
                            stat[17]++;
                            break;
                        case 'D':
                            stat[18]++;
                            break;
                        case 'E':
                            stat[19]++;
                            break;
                        }
                    }

                    i = -2; //not very good idea
                    j = 0;
                    seq_num++;
                }
                i++;
            }

            total = 120 * seq_num;

            for (k = 0; k < 20; k++) {
                distribution[k] = (float) 1000.0 * stat[k] / total;
                switch (k) {
                case 0:
                    System.out.print("Amino Acid G: ");
                    System.out.print(distribution[0]);
                    System.out.print("\n");
                    break;
                case 1:
                    System.out.print("Amino Acid A: ");
                    System.out.print(distribution[1]);
                    System.out.print("\n");
                    break;
                case 2:
                    System.out.print("Amino Acid P: ");
                    System.out.print(distribution[2]);
                    System.out.print("\n");
                    break;
                case 3:
                    System.out.print("Amino Acid V: ");
                    System.out.print(distribution[3]);
                    System.out.print("\n");
                    break;
                case 4:
                    System.out.print("Amino Acid L: ");
                    System.out.print(distribution[4]);
                    System.out.print("\n");
                    break;
                case 5:
                    System.out.print("Amino Acid I: ");
                    System.out.print(distribution[5]);
                    System.out.print("\n");
                    break;
                case 6:
                    System.out.print("Amino Acid M: ");
                    System.out.print(distribution[6]);
                    System.out.print("\n");
                    break;
                case 7:
                    System.out.print("Amino Acid F: ");
                    System.out.print(distribution[7]);
                    System.out.print("\n");
                    break;
                case 8:
                    System.out.print("Amino Acid Y: ");
                    System.out.print(distribution[8]);
                    System.out.print("\n");
                    break;
                case 9:
                    System.out.print("Amino Acid W: ");
                    System.out.print(distribution[9]);
                    System.out.print("\n");
                    break;
                case 10:
                    System.out.print("Amino Acid S: ");
                    System.out.print(distribution[10]);
                    System.out.print("\n");
                    break;
                case 11:
                    System.out.print("Amino Acid T: ");
                    System.out.print(distribution[11]);
                    System.out.print("\n");
                    break;
                case 12:
                    System.out.print("Amino Acid C: ");
                    System.out.print(distribution[12]);
                    System.out.print("\n");
                    break;
                case 13:
                    System.out.print("Amino Acid N: ");
                    System.out.print(distribution[13]);
                    System.out.print("\n");
                    break;
                case 14:
                    System.out.print("Amino Acid Q: ");
                    System.out.print(distribution[14]);
                    System.out.print("\n");
                    break;
                case 15:
                    System.out.print("Amino Acid K: ");
                    System.out.print(distribution[15]);
                    System.out.print("\n");
                    break;
                case 16:
                    System.out.print("Amino Acid H: ");
                    System.out.print(distribution[16]);
                    System.out.print("\n");
                    break;
                case 17:
                    System.out.print("Amino Acid R: ");
                    System.out.print(distribution[17]);
                    System.out.print("\n");
                    break;
                case 18:
                    System.out.print("Amino Acid D: ");
                    System.out.print(distribution[18]);
                    System.out.print("\n");
                    break;
                case 19:
                    System.out.print("Amino Acid E: ");
                    System.out.print(distribution[19]);
                    System.out.print("\n");
                    break;
                }
            }

            infile.close();
        } catch (IOException e) {
           throw e;
        }
    }
}
