/**
  \file        doc_rotator.cpp
  \author      Arnaud Ramey <arnaud.a.ramey@gmail.com>
                -- https://sites.google.com/site/rameyarnaud
  \date        2012/07

________________________________________________________________________________

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
________________________________________________________________________________

A simple tool for rectifying tilted scanned images.

It first detects the lines segments in the document.
These lines can be made by text characters, drawings, etc.
The method used is the Hough transform for line segments, cf the Wikipedia page about it.
http://en.wikipedia.org/wiki/Hough_transform

However, not all lines correspond to the text lines, there also are some
outliers, such as a frame around the page, or simply false detections. To
determine the predominant orientation, a RANSAC filter is applied on the
bunch of line segments (cf the RANSAC Wikipedia page for details:
http://en.wikipedia.org/wiki/RANSAC ).

The most "popular" orientation is then used to rectify the document.
 */


#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <stdio.h>

#include <ransac.h>

inline double absolute_angle_between_two_angles(const double & a1,
                                                const double & a2) {
  double lowest_error = 10;
  for (int mod = -1; mod <= 1; ++mod) {
    double curr_error = fabs(a1 + mod * 2 * PI - a2);
    if (curr_error < lowest_error)
      lowest_error = curr_error;
  } // end loop mod
  return lowest_error;
}

////////////////////////////////////////////////////////////////////////////////

/*!
 * \brief   returns the average angle of a bunch of points,
            taking care of the painful case around [+-pi]
 * \param   angles some angles in the interval [-pi,+pi] radians
 * \return  the average angle , in the interval [-pi,+pi] radians.
 */
static inline double oriented_average_angle(const std::vector<double> & angles,
                                            double & avg_dist) {
  // compute average
  double x = 0, y = 0;
  for (unsigned int angle_idx = 0; angle_idx < angles.size(); ++angle_idx) {
    x += cos(angles[angle_idx]);
    y += sin(angles[angle_idx]);
  } // end loop angle_idx
  x /= angles.size();
  y /= angles.size();
  double avg_angle = atan2(y, x);

  // determine distance
  avg_dist = 0;
  // examinate angle - 2 PI, angle, angle + 2 PI
  for (unsigned int angle_idx = 0; angle_idx < angles.size(); ++angle_idx)
    avg_dist += absolute_angle_between_two_angles(angles[angle_idx], avg_angle);
  avg_dist /= angles.size();

  return avg_angle;
}

////////////////////////////////////////////////////////////////////////////////

static inline void rotate_img_from_pt(const cv::Mat & src,
                                      cv::Mat & dst,
                                      const cv::Point & center,
                                      const double angle_rad,
                                      const double scale = 1.,
                                      cv::Scalar border_value = cv::Scalar::all(255)) {
  cv::Mat rot_mat = getRotationMatrix2D( center, angle_rad * 180.f / M_PI, scale );
  cv::warpAffine( src, dst, rot_mat, src.size(),
                  cv::INTER_LINEAR, cv::BORDER_CONSTANT,
                  border_value);
}

////////////////////////////////////////////////////////////////////////////////

class MeanModel : public ransac::Model<double>{
public:

  //! see Model::fit_to_vector()
  double fit_to_vector(const DataSet & data) {
    //    printf("\nfit_to_vector(%s)",
    //                 StringUtils::iterable_to_string(data).c_str());

    double avg_dist;
    average_angle = oriented_average_angle(data, avg_dist);
    return avg_dist;
  } // end fit_to_vector()

  //! see Model::distance_to_current()
  double distance_to_current(const double & query) const {
    return absolute_angle_between_two_angles(average_angle, query);
  } // end distance_to_current()

  //private:
  double average_angle;
};

////////////////////////////////////////////////////////////////////////////////

inline cv::Mat process_file(const cv::Mat & img_in) {
  //printf("process_file()\n");
  // color -> b&w
  cv::Mat1b img_edges;
  cv::cvtColor(img_in, img_edges, cv::COLOR_BGR2GRAY);
  cv::Canny( img_edges, img_edges, 50, 200, 3 );
  // resize
  double factor = 1000.f / img_edges.cols;
  cv::resize(img_edges, img_edges, cv::Size(), factor, factor);
  // hough
  std::vector<cv::Vec4i> lines;
  cv::HoughLinesP(img_edges,
                  lines,
                  1, //Distance resolution in pixel-related units
                  CV_PI / 180.f, //Angle resolution measured in radians
                  300, // Accumulator threshold parameter. Only those lines are returned that get enough votes (  ).
                  100, // Minimum line length. Line segments shorter than that are rejected.
                  100); // Maximum allowed gap between points on the same line to link them.
  printf("\nThere are %i lines.\n", lines.size());

  // draw lines
  //cv::Mat3b img_illus = img_in.clone();
  cv::Mat3b img_illus;
  cv::cvtColor(img_edges, img_illus, cv::COLOR_GRAY2BGR);
  for (unsigned int line_idx = 0; line_idx < lines.size(); ++line_idx) {
    // draw blue lines
    cv::line(img_illus,
             cv::Point(lines[line_idx][0], lines[line_idx][1]),
             cv::Point(lines[line_idx][2], lines[line_idx][3]),
             cv::Scalar(255, 0, 0), 1);
  } // end loop line_idx

  // compute angles
  std::vector<double> angles;
  for (unsigned int line_idx = 0; line_idx < lines.size(); ++line_idx) {
    angles.push_back(atan2(
                       lines[line_idx][3] - lines[line_idx][1],
                       lines[line_idx][2] - lines[line_idx][0]
                       ));
    // printf("%f\t", angles.back() * 180.f / M_PI);
  }
  // printf("\n");

  // compute average angle with ransac
  MeanModel initial, answer;
  double best_error = 0;
  std::vector<double> best_consensus;
  ransac::ransac<double, MeanModel>(angles, initial, 5, 50, .1, 10,
                                    answer, best_consensus, best_error);
  double consensus_average_angle_rad = answer.average_angle;
  double consensus_average_angle_deg = consensus_average_angle_rad * 180.f / M_PI;
  printf("best_consensus: size %i, best_error:%f, "
         "consensus_average_angle_rad:%f, "
         "consensus_average_angle_deg:%f\n",
         best_consensus.size(), best_error,
         consensus_average_angle_rad, consensus_average_angle_deg);
  if (fabs(sin(consensus_average_angle_rad)) > 0.2 ) {
    printf("angle too remote from 0, doing nothing\n");
    return img_in;
  }

  // print best_consensus
  //  for (unsigned int cons_idx = 0; cons_idx < best_consensus.size(); ++cons_idx)
  //    printf("%g\t", best_consensus[cons_idx]);
  //  printf("\n");

  // draw compass on img_illus
  double compass_size = img_illus.cols / 2;
  cv::Point pt1(compass_size, compass_size);
  cv::Point pt2(pt1.x + cos(consensus_average_angle_rad) * compass_size,
                pt1.y + sin(consensus_average_angle_rad) * compass_size);
  cv::line(img_illus, pt1, pt2, cv::Scalar(0, 0, 255), 3);

  // rotate image
  cv::Mat img_out;
  rotate_img_from_pt(img_in, img_out,
                     cv::Point(img_in.cols / 2, img_in.rows / 2),
                     consensus_average_angle_rad);


  // show images
#if 0 // some debug stuff
  // scale images down
  //cv::imshow("img_edges", img_edges);
  cv::Mat3b img_in_copy;
  cv::resize(img_in, img_in_copy, cv::Size(300, img_in.rows * 300 / img_in.cols));
  cv::resize(img_illus, img_illus, cv::Size(300, img_illus.rows * 300 / img_illus.cols));
  cv::resize(img_out, img_out, cv::Size(300, img_out.rows * 300 / img_out.cols));
  cv::imshow("img_in", img_in_copy);
  cv::imshow("img_illus", img_illus);
  cv::imshow("img_out", img_out);
  cv::waitKey(0);
#else
  cv::imshow("img_in", img_in);
  cv::imshow("img_illus", img_illus);
  cv::imshow("img_out", img_out);
  cv::waitKey(25);
#endif
  return img_out;
} // end process_file();

////////////////////////////////////////////////////////////////////////////////


/*!
  Add a suffix to a filename
 \param path
    The absolute or relative path to filename
 \param suffix
    The string to be added to path, before the file extension
 \return std::string
 \example ("/foo/bar.dat", "_out") returns "/foo/bar_out.dat"
          ("/foo/bar", "_out") returns "/foo/bar_out"
*/
inline std::string add_suffix_before_filename_extension
(const std::string & path, const std::string & suffix = "out") {
  std::string::size_type dot_pos = path.find_last_of('.');
  std::ostringstream out;
  if (dot_pos == std::string::npos) {
    out << path << suffix;
  }
  else {
    std::string path_before_dot = path.substr(0, dot_pos);
    std::string path_after_dot = path.substr(dot_pos + 1);
    out << path_before_dot << suffix << "." << path_after_dot;
  }
  return out.str();
}

////////////////////////////////////////////////////////////////////////////////

inline void process_file(std::string filename) {
  printf("\nprocess_file('%s')\n", filename.c_str());
  cv::Mat in = cv::imread(filename, CV_LOAD_IMAGE_COLOR);
#if 0
  // crook image for tests
  cv::Mat in2;
  rotate_img_from_pt(in, in2, cv::Point(in.cols / 2, in.rows / 2), (2 * drand48() - 1) * .1);
  in2.copyTo(in);
#endif

  cv::Mat out = process_file(in);
  std::string filename_out = add_suffix_before_filename_extension(filename, "_out");
  cv::imwrite(filename_out, out);
  printf("Image saved in %s'\n", filename_out.c_str());
} // end process_file();

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
  printf("\ndoc_rotator\n");
  srand(time(NULL));
  if (argc < 2) {
    printf("\nuse: '%s' [filename]\n", argv[0]);
    return -1;
  }
  for (unsigned int file_idx = 1; file_idx < argc; ++file_idx) {
    std::string filename = argv[file_idx];
    process_file(filename);
  } // end loop file_idx
  return 0;
}
