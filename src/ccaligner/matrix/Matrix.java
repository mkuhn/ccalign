/*
 * $Id: Matrix.java,v 1.2 2005/04/14 14:44:42 ahmed Exp $
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

package ccaligner.matrix;

import java.io.Serializable;

/**
 * Scoring matrix.
 * 
 * @author Ahmed Moustafa (ahmed@users.sf.net)
 */

public class Matrix implements Serializable {
    /**
     * 
     */
    private static final long serialVersionUID = 3833742170619524400L;

    /**
     * Matrix id (or name)
     */
    private String id = null;
    
    /**
     * Scores
     */
    private float[][] scores = null;
    
    
    public Matrix() {
        this.id = "zeroes";
        this.scores = new float[127][127];
    }
    
    public Matrix(String id, float[][] scores) {
        this.id = id;
        this.scores = scores;
    }
    
    /**
     * @return Returns the id.
     */
    public String getId() {
        return this.id;
    }
    
    /**
     * @return Check scores.
     */
    public boolean checkScores() {
    	char[] aa = "ARNDCEQGHILKMFPSTWYV".toCharArray();
    	for (char c1 : aa)
    	{
    		float max = 0;
        	for (char c2 : aa)
    		{
        		float f = this.scores[c1][c2];
    			if (f > max) { max = f; };
    		}
        	if (id.contentEquals("zeroes"))
        	{
        		assert max == 0 : this.id + ": " + c1;
        		if (max > 0) return false; 
        	}
        	else
        	{
	    		assert max > 0 : this.id + ": " + c1;
	    		if (max <= 0) return false;
	        }
    	}
        return true;
    }

    /**
     * @return Returns the scores.
     */
    public float[][] getScores() {
    	assert checkScores();
        return scores;
    }
    
    public void scaleScores(float f)
    {
    	for (int i=0; i<127; i++)
    	{
        	for (int j=0; j<127; j++)
        	{
        		scores[i][j] *= f;
        	}
    	}
    }

    public void roundScores()
    {
    	for (int i=0; i<127; i++)
    	{
        	for (int j=0; j<127; j++)
        	{
        		scores[i][j] = Math.round(scores[i][j]);
        	}
    	}
    }

    /**
     * 
     * @param a
     * @param b
     * @return score
     */
    public float getScore(char a, char b) {
        return this.scores[a][b];
    }
}