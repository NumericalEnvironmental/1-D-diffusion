# 1-D-diffusion
This script, written in julia, simulates diffusion of a gaseous constituent through a partially water-saturated 1-D column with fixed concentration boundaries at either end (one is set equal to zero). An analytical solution for 1-D diffusion into a semi-infinite half-space is also included for comparison. This script will later be incorporated into a more flexible geometry reactive transport model as par tof an ongoing project.

Two text input files are included:

parameters.txt --> characteristics of the 1-D model domain as well as the gaseous species

knobs.txt --> miscellaneous constraints in the numerical solver

More info can be found here: https://numericalenvironmental.wordpress.com/2016/11/23/two-phase-reactive-transport-model-development-with-julia-part-1-a-simple-1-d-diffusion-model/

I'd appreciate hearing back from you if you find the script useful. Questions or comments are welcome at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
