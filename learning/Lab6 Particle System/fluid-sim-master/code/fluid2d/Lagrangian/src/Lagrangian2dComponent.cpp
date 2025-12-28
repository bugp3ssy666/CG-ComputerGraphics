#include "Lagrangian2dComponent.h"

#include "Lagrangian/include/Lagrangian2dComponent.h"

namespace FluidSimulation
{
    namespace Lagrangian2d
    {
        // ========== 喷泉场景参数 ==========
        // 喷泉源位置（归一化坐标，容器中心底部）
        static glm::vec2 fountainSourcePos = glm::vec2(0.0f, -0.8f);
        // 喷泉初始速度（向上喷出）
        static glm::vec2 fountainVelocity = glm::vec2(0.0f, 2.5f);
        // 每次模拟步骤生成的粒子数量
        static int particlesPerStep = 3;
        // 喷泉水平扩散范围
        static float fountainSpread = 0.1f;
        // 是否启用喷泉模式
        static bool enableFountain = true;
        // 当前帧计数器（用于控制生成频率）
        static int frameCounter = 0;

        void Lagrangian2dComponent::shutDown()
        {
            delete renderer, solver, ps;
            renderer = NULL;
            solver = NULL;
            ps = NULL;
            frameCounter = 0;
        }

        void Lagrangian2dComponent::init()
        {
            if (renderer != NULL || solver != NULL || ps != NULL)
            {
                shutDown();
            }

            Glb::Timer::getInstance().clear();

            // ========== 初始化渲染器 ==========
            renderer = new Renderer();
            renderer->init();

            // ========== 初始化粒子系统 ==========
            // 设置容器大小
            ps = new ParticleSystem2d();
            ps->setContainerSize(glm::vec2(-1.0f, -1.0f), glm::vec2(1.0f, 1.0f));

            // ========== 设置喷泉场景 ==========
            // 如果配置中有流体块，添加它们
            // 否则，喷泉模式将在simulate中持续添加粒子
            for (int i = 0; i < Lagrangian2dPara::fluidBlocks.size(); i++) {
                ps->addFluidBlock(Lagrangian2dPara::fluidBlocks[i].lowerCorner, Lagrangian2dPara::fluidBlocks[i].upperCorner, 
                    Lagrangian2dPara::fluidBlocks[i].initVel, Lagrangian2dPara::fluidBlocks[i].particleSpace);
            }

            // 如果没有预设的流体块，初始生成一些喷泉粒子
            if (Lagrangian2dPara::fluidBlocks.size() == 0 && enableFountain)
            {
                // 初始喷出一批粒子
                ps->addFountainParticles(fountainSourcePos, fountainVelocity, 50, fountainSpread);
            }

            ps->updateBlockInfo();

            Glb::Logger::getInstance().addLog("2d Particle System initialized (Fountain Mode). Initial particle num: " + std::to_string(ps->particles.size()));

            solver = new Solver(*ps);
            frameCounter = 0;
        }

        void Lagrangian2dComponent::simulate()
        {
            // ========== 喷泉粒子管理 ==========
            // 每隔几帧生成新粒子，控制喷泉流量
            frameCounter++;
            if (enableFountain && frameCounter % 2 == 0)  // 每2帧生成一次
            {
                // 在喷泉源位置生成新粒子
                ps->addFountainParticles(fountainSourcePos, fountainVelocity, particlesPerStep, fountainSpread);
            }

            // 移除超出边界的死亡粒子（可选，保持粒子数量在合理范围）
            if (frameCounter % 50 == 0)  // 每50帧清理一次
            {
                ps->removeDeadParticles();
            }

            // ========== SPH求解循环 ==========
            for (int i = 0; i < Lagrangian2dPara::substep; i++)
            {
                ps->updateBlockInfo();
                solver->solve();
            }
        }

        GLuint Lagrangian2dComponent::getRenderedTexture()
        {
            renderer->draw(*ps);
            return renderer->getRenderedTexture();
        }
    }
}
