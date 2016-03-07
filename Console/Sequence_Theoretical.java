package Console;

public class Sequence_Theoretical {
	public static void main(String[] arg) {
		int aa[] = new int[20];
		aa[0] = 72;
		aa[1] = 78;
		aa[2] = 52;
		aa[3] = 66;
		aa[4] = 91;
		aa[5] = 53;
		aa[6] = 23;
		aa[7] = 39;
		aa[8] = 32;
		aa[9] = 14;
		aa[10] = 68;
		aa[11] = 59;
		aa[12] = 19;
		aa[13] = 43;
		aa[14] = 42;
		aa[15] = 59;
		aa[16] = 23;
		aa[17] = 51;
		aa[18] = 53;
		aa[19] = 63;
		int value = 0;
		int i = 0, j = 0;
		for (i = 0; i < 20; i++) {
			for (j = 0; j < 20; j++) {
				value = aa[j] * aa[i] * 119;
				System.out.print(value + "\n");
			}
		}
	}
}
