package Visualize;

import java.awt.image.*;

import javax.imageio.*;

import java.awt.*;
import java.io.File;
import java.io.IOException;

public class Pic_test {

	public static void main(String[] args) throws IOException {
		BufferedImage bufferedImage = new BufferedImage(1026, 1026, BufferedImage.TYPE_INT_RGB);
	    
        // Create a graphics contents on the buffered image
        Graphics2D g2d = bufferedImage.createGraphics();
        Color co = new Color(255, 0, 0);
		g2d.setColor(co);
		g2d.fillRect(1023, 1024, 1024, 1024);
		try {
            // Save as PNG
            File file = new File("D:\\Pic_try.png");
            ImageIO.write(bufferedImage, "png", file);
        } catch (IOException e) {
        	throw e;
        }

	}

}
