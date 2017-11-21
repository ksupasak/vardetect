/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotec.bsi.ngs.vardetect.core.util;

import java.net.URL;
import java.net.URLConnection;
import java.net.HttpURLConnection;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.io.Reader;

/**
 *
 * @author worawich
 */
public class EnsemblAPI {
    public static String getSequence(String species, String chromosomeName, long startPosition, long stopPosition) throws Exception {
        /**
         * This function use to get portion of sequence from ensembl online database by specific species, chromosome, start and stop position
         */
        
        String server = "http://rest.ensembl.org";
        String start = Long.toString(startPosition);
        String stop = Long.toString(stopPosition);
        String ext = "/sequence/region/"+species+"/"+chromosomeName+":"+startPosition+".."+stopPosition+":1?";
        URL url = new URL(server + ext);

        URLConnection connection = url.openConnection();
        HttpURLConnection httpConnection = (HttpURLConnection)connection;

        httpConnection.setRequestProperty("Content-Type", "text/plain");

        InputStream response = connection.getInputStream();
        int responseCode = httpConnection.getResponseCode();

        if(responseCode != 200) {
          throw new RuntimeException("Response code was not 200. Detected response was "+responseCode);
        }
 
        String output;
        Reader reader = null;
        try {
          reader = new BufferedReader(new InputStreamReader(response, "UTF-8"));
          StringBuilder builder = new StringBuilder();
          char[] buffer = new char[8192];
          int read;
          while ((read = reader.read(buffer, 0, buffer.length)) > 0) {
            builder.append(buffer, 0, read);
          }
          output = builder.toString();
        } 
        finally {
            if (reader != null) try {
              reader.close(); 
            } catch (IOException logOrIgnore) {
              logOrIgnore.printStackTrace();
            }
        }
        return output;
    } 
}
