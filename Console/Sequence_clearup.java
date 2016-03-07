package Console;

/*Step 1
 *remove incomplete sequences
 *remove poly G (caused by server connection failure), five or more continues G is not allowed
 *file obtained: std_aa_1M_queried.xml, contains 1000000 sequences
 *use this program carefully, as will take a lot of time
 */
import java.io.*;

public class Sequence_clearup {
    public static void main(String[] args) throws IOException {

        String infilename =
                "D:\\My Documents\\GThesis\\DatabaseBAK\\std_aa_1M_02.xml";
        String outfilename =
                "D:\\My Documents\\GThesis\\DatabaseBAK\\std_aa_1M_02_queried.xml";

        byte data[] = new byte[1];
        byte line[] = new byte[140];
        int i = 0, G = 0;
        int error = 0;
        int polyG_error = 0;
        long useful = 0;

        try {
            FileInputStream infile = new FileInputStream(infilename);
            FileOutputStream outfile = new FileOutputStream(outfilename);

            while (((infile.read(data, 0, 1)) != -1) && (useful < 1000000)) {
                if (data[0] != 10){
                    line[i] = data[0];
                    i ++;
                }else if(i == 139){
                	//remove poly G
                	for (G = 0; G < i; G++) {
                		if (line[G] == 71) {
                			G ++;
                			if (line[G] == 71) {
                				G ++;
                    			if (line[G] == 71) {
                    				G ++;
                        			if (line[G] == 71) {
                        				G ++;
                                		if (line[G] == 71) {
                                			polyG_error ++;
                                			System.out.print("polyG error detected!\n");
                                			break;
                                		}else
                                			continue;
                        			}else
                        				continue;
                    			}else
                    				continue;
                			}else
                				continue;
                		}else
                			continue;
                	}
                	if (G == i) {	                	
	                    System.out.print("sequence passed, ready to input \n");
	                    outfile.write(line, 0, 138);
	                    useful ++;
	                    outfile.write(13);  //carriage return
	                    outfile.write(10);  //new line
                	}
                    i = 0;
                }else{
                    System.out.print("error \n");
                    error ++;
                    i = 0;
                }
            }
            System.out.print("Done\n");
            System.out.print("error found : " + error + "\n");
            System.out.print("PolyG error found : " + polyG_error + "\n");
            System.out.print("Useful Sequence Number : " + useful + "\n");
        } catch (IOException e) {
            throw e;
        }
    }
}
